#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <optional>
#include "nanojpeg.hpp"
#include "mio.hpp"
#include <algorithm>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "turbojpeg.h"

#ifndef stbi__float2fixed
#define stbi__float2fixed(x) (((int)((x) * 4096.0f + 0.5f)) << 8)
#endif
template <typename T>
constexpr inline void YCbCr_to_RGB(int y, int cb, int cr, T &r_, T &g_, T &b_)
{
    int y_fixed = (y << 20) + (1 << 19); // rounding
    int r{}, g{}, b{};
    cr -= 128;
    cb -= 128;
    r = y_fixed + cr * stbi__float2fixed(1.40200f);
    g = y_fixed + (cr * -stbi__float2fixed(0.71414f)) + ((cb * -stbi__float2fixed(0.34414f)) & 0xffff0000);
    b = y_fixed + cb * stbi__float2fixed(1.77200f);
    r >>= 20;
    g >>= 20;
    b >>= 20;
    if ((unsigned)r > 255)
    {
        if (r < 0)
            r = 0;
        else
            r = 255;
    }
    if ((unsigned)g > 255)
    {
        if (g < 0)
            g = 0;
        else
            g = 255;
    }
    if ((unsigned)b > 255)
    {
        if (b < 0)
            b = 0;
        else
            b = 255;
    }
    r_ = r;
    g_ = g;
    b_ = b;
}

template <bool is_ycck, int planes_nb, typename YUV_, typename RGB_, typename COMP>
inline void convert(YUV_ &&yuv, RGB_ &&rgb, COMP &&comp, int width, int height)
{

    for (int h = 0; h < height; ++h)
        for (int x = 0; x < width; ++x)
        {
            if constexpr (planes_nb == 1)
            {
                rgb(x, h, std::clamp(yuv(0, x, h), 0, 255));
            }
            else if constexpr (planes_nb == 2)
            {
                rgb(x, h, std::clamp(yuv(0, x, h), 0, 255), std::clamp(yuv(1, x, h), 0, 255));
            }
            else if constexpr (planes_nb == 3)
            {
                auto y = yuv(0, x, h);
                auto u = yuv(1, x >> comp[1].chroma_w_log2, h >> comp[1].chroma_h_log2);
                auto v = yuv(2, x >> comp[2].chroma_w_log2, h >> comp[2].chroma_h_log2);
                int r, g, b;
                YCbCr_to_RGB(y, u, v, r, g, b);
                rgb(x, h, r, g, b);
            }
            else
            {
                const auto c = yuv(0, x >> comp[0].chroma_w_log2, h >> comp[0].chroma_h_log2);
                const auto m = yuv(1, x >> comp[1].chroma_w_log2, h >> comp[1].chroma_h_log2);
                const auto y = yuv(2, x >> comp[2].chroma_w_log2, h >> comp[2].chroma_h_log2);
                const auto k = yuv(3, x >> comp[3].chroma_w_log2, h >> comp[3].chroma_h_log2);
                if constexpr (is_ycck)
                {
                    int ir, ig, ib;
                    YCbCr_to_RGB(c, m, y, ir, ig, ib);
                    const auto r = nanojpeg::njClip((255 - ir) * k / 255);
                    const auto g = nanojpeg::njClip((255 - ig) * k / 255);
                    const auto b = nanojpeg::njClip((255 - ib) * k / 255);
                    rgb(x, h, r, g, b);
                }
                else
                {
                    const auto r = nanojpeg::njClip(c * k / 255);
                    const auto g = nanojpeg::njClip(m * k / 255);
                    const auto b = nanojpeg::njClip(y * k / 255);
                    rgb(x, h, r, g, b);
                }
            }
        }
}

template <typename YUV_, typename RGB_>
inline void convert(nanojpeg::nj_result &image, YUV_ &&yuv, RGB_ &&rgb)
{
    switch (image.components.size())
    {
    case 1:
        convert<false, 1>(std::forward<YUV_>(yuv), std::forward<RGB_>(rgb), image.components, image.width, image.height);
        break;
    case 2:
        convert<false, 2>(std::forward<YUV_>(yuv), std::forward<RGB_>(rgb), image.components, image.width, image.height);
        break;
    case 3:
        convert<false, 3>(std::forward<YUV_>(yuv), std::forward<RGB_>(rgb), image.components, image.width, image.height);
        break;
    case 4:
        if (image.is_ycck)
            convert<true, 4>(std::forward<YUV_>(yuv), std::forward<RGB_>(rgb), image.components, image.width, image.height);
        else
            convert<false, 4>(std::forward<YUV_>(yuv), std::forward<RGB_>(rgb), image.components, image.width, image.height);
        break;
    default:
        throw std::runtime_error("invalid number of planes");
    }
}

template <typename _STR, typename COMPS>
void save_bmp(_STR &&out_filename, COMPS &&image)
{

    int comp_nb = std::min<int>(3, image.components.size());
    std::vector<uint8_t> rgb(image.width * image.height * comp_nb);
    auto rgb_out = rgb.data();
    convert(image, [&comp = std::forward<COMPS>(image.components)](int comp_n, int x, int y)
            { return (int)comp[comp_n].pixels[y * comp[comp_n].stride + x]; }, [&rgb_out](int x, int y, auto &&...args)
            { ((*rgb_out++ = args), ...); });
    stbi_write_bmp(out_filename.c_str(), image.width, image.height, comp_nb, rgb.data());
}

namespace profiling
{
    using std::chrono::duration;
    using std::chrono::high_resolution_clock;

    struct StopWatch
    {
        std::optional<high_resolution_clock::time_point> start_point{};
        double elapsed_total{};

        inline void startnew()
        {
            elapsed_total = 0.0;
            start();
        }

        inline void start()
        {
            if (!start_point)
                start_point = high_resolution_clock::now();
        }

        inline bool is_running() const { return start_point.has_value(); }

        inline double elapsed() const
        {
            return elapsed_total + ((start_point) ? duration<double>(high_resolution_clock::now() - *start_point).count() : 0);
        }
        void stop()
        {
            if (start_point)
            {
                elapsed_total += duration<double>(high_resolution_clock::now() - *start_point).count();
                start_point.reset();
            }
        }
        std::string elapsed_str() const
        {
            return std::string((std::stringstream() << std::fixed << std::setprecision(9) << elapsed()).str());
        };
    };
}
using profiling::StopWatch;

struct comp
{
    std::vector<uint8_t> data{};
    int chroma_h_log2{}, chroma_w_log2{};
    void setup(int subsamp, int comp)
    {
        switch (comp)
        {
        case 0:
            chroma_h_log2 = 0;
            chroma_w_log2 = 0;
            break;
        case 1:
        case 2:
            switch (subsamp)
            {
            case TJSAMP_444:
                chroma_h_log2 = 0;
                chroma_w_log2 = 0;
                break;
            case TJSAMP_420:
                chroma_h_log2 = 1;
                chroma_w_log2 = 1;
                break;
            case TJSAMP_422:
                chroma_h_log2 = 0;
                chroma_w_log2 = 1;
                break;
            }
        }
    }
};

__attribute__((noinline, optnone))
static auto jpeg_turbo_bench(StopWatch &bench, const uint8_t *buf, size_t size, bool fastDCT)
{

    nanojpeg::nj_context_t njheader{};
    njheader.pos = buf;
    njheader.size = size;
    njheader.DecodeSOF<true, false>();
    if (njheader.ncomp != 3)
        throw std::runtime_error("jpeg turbo can only decode direct YUV8 places");


    struct H
    {
        tjhandle handle{};
        H(tjhandle h) : handle(h) {}
        ~H() { tj3Destroy(handle); }
    } tj3handleRAII(tj3Init(TJINIT_DECOMPRESS));

    // tj3DecompressHeader(tj3handleRAII.handle, buf, size); //without this call we can not retrive the buffer size
    njheader.allocate_pixels(); // bench memory on tj side

    uint8_t *yuv_planes[3] = {njheader.comp[0].pixels.data(), njheader.comp[1].pixels.data(), njheader.comp[2].pixels.data()};
    int strides[3] = {njheader.comp[0].stride, njheader.comp[1].stride, njheader.comp[2].stride};

    if (fastDCT)
        tj3Set(tj3handleRAII.handle, TJPARAM_FASTDCT, 1);
    bench.start();
    if (tj3DecompressToYUVPlanes8(tj3handleRAII.handle, buf, size, (uint8_t **)yuv_planes, (int *)strides) != 0)
        throw std::runtime_error("tjDecompressToYUVPlanes failed:" + std::string(tjGetErrorStr2(tj3handleRAII.handle)));
    bench.stop();

    return njheader.comp;
}

__attribute__((noinline, optnone))
static auto nanojpeg_bench(StopWatch &bench, const uint8_t *buf, size_t size)
{
    nanojpeg::nj_context_t njheader{};
    njheader.pos = buf;
    njheader.size = size;
    njheader.DecodeSOF<true, false>();
    nanojpeg::nj_result frame{};
    //avoide memory allocation in nanojpeg side
    frame.components.resize( njheader.ncomp );
    for (int i = 0; i < njheader.ncomp; i++) {
        frame.components[i].pixels.resize( njheader.comp[i].size);
    }

    bench.start();
    nanojpeg::decode(buf, size, frame);
    bench.stop();
    return frame;
}

__attribute__((noinline, optnone))
static void nanojpeg_motion_bench(StopWatch &bench, const uint8_t *&buf, size_t &size, nanojpeg::nj_result &frame)
{
    bench.start();
    nanojpeg::decode(buf, size, frame);
    bench.stop();
}


inline void print_stat(std::string title, double elapsed, double imgMPixSize, size_t size_total, int times)
{
    std::cout << "** " << title << " **" << std::endl;
    std::cout << "time   = " << std::fixed << std::setprecision(9) << (elapsed / times) << " seconds" << std::endl;
    std::cout << "images = " << std::fixed << std::setprecision(3) << (times / elapsed) << " fps" << std::endl;
    std::cout << "bytes  = " << std::fixed << std::setprecision(2) << (size_total / (1024.0 * 1024.0) / elapsed) << " MB/s" << std::endl;
    std::cout << "pixels = " << std::fixed << std::setprecision(2) << (imgMPixSize * double(times) / elapsed) << " MPix/s" << std::endl;
    std::cout << std::endl;
}

template <typename convertable_to_string>
inline std::string ext(convertable_to_string &&path)
{
    auto s = std::string(path);
    auto pos = s.find_last_of('.');
    if ( std::string::npos == pos) {
        return "";
    }
    else {
        std::string ext = s.substr(pos + 1);
        for (auto &c : ext)
            c = std::tolower(c);
        return ext;
    }
}

int main(int argc, char **argv)
{

    try
    {
        std::error_code error;
        mio::ummap_source mmap = mio::make_mmap<mio::ummap_source>(argv[1], 0, 0, error);
        if (error)
        {
            throw std::runtime_error(error.message());
        }
        StopWatch njtime{}, tjtime{}, tjtime_fast{}, convert_time{};
        auto file_extenstion = ext(argv[1]);
        bool is_motion = file_extenstion == "mjpeg" ||  file_extenstion == "mjpg";
        if ( is_motion ) {
            std::cout << "Motion JPEG detected" << std::endl;
            StopWatch motion_time{};

            const uint8_t * pos = mmap.data();
            size_t    size = mmap.size();
            int times = 0; // number of frames decoded
            std::ofstream fileyuv_file( std::string(argv[1]) + ".y4m", std::ios::binary | std::ios::out | std::ios::trunc);
            nanojpeg::nj_result frame{};

            while (size > 0)
            {
                motion_time.start();
                nanojpeg_motion_bench(motion_time, pos, size,frame);
                motion_time.stop();
                times+=1;
                if (times == 1) {
                    //header
                    if ( frame.yuv_format != 420 )
                        fileyuv_file << "YUV4MPEG2 W" << frame.width << " H" << frame.height << " F25:1 It C" << frame.yuv_format <<  " XYSCSS=" << frame.yuv_format << " XCOLORRANGE=FULL\n";
                    else
                        fileyuv_file << "YUV4MPEG2 W" << frame.width << " H" << frame.height << " F25:1 It C420jpeg XYSCSS=420JPEG XCOLORRANGE=FULL\n";

                }
                if (times % 100 == 0)
                    std::cout << "." << std::flush;

                fileyuv_file << "FRAME\n";

                for (const auto& c : frame.components)
                {
                    const uint8_t * p = c.pixels.data();
                    for (int h = 0; h < frame.height >> c.chroma_h_log2; h++, p+=c.stride)
                    {
                        fileyuv_file.write(reinterpret_cast<const char *>(p), frame.width >> c.chroma_w_log2 );
                    }
                }
                fileyuv_file.flush();
            }
            std::cout << std::endl;
            auto imgMPixSize = frame.width * frame.height * 1e-6;
            std::cout << "stream  = " << frame.width << "x" << frame.height << " (" << std::fixed << std::setprecision(1) << imgMPixSize << " MPix)" << std::endl;

            print_stat("nanojpeg motion", motion_time.elapsed(), imgMPixSize,  mmap.size(), times);
            return 0;
        }

        auto image = nanojpeg::decode(mmap.data(), mmap.size());//load memmap into nanojpeg image
        if (image.yuv_format == 0) {
            std::cerr << "warining: not standard yuv format!" << std::endl;
        } else {
            std::cout << "* YUV" << image.yuv_format << std::endl;
            std::cout << "* " << image.width << "x" << image.height << std::endl;
        }

        int times = std::clamp(1024 * 1024 * 1024.0 / image.size,1.,2000.);
        std::cout << "Benchmarking repeats " << times << " times. Every dot is 100 frames." << std::endl;

        bool turbo_ok = true;
        try
        {
            StopWatch dummy{};
            jpeg_turbo_bench(dummy, mmap.data(), mmap.size(), true);
        }
        catch (const std::exception &e)
        {
            std::cout << std::endl
                      << "turbojpeg failed: " << e.what() << std::endl;
            turbo_ok = false;
        }

        for (int i = 0; i < times; i++)
        {

            (void)nanojpeg_bench(njtime, mmap.data(), mmap.size());

            if (turbo_ok)
            {
                (void)jpeg_turbo_bench(tjtime, mmap.data(), mmap.size(), false);
                (void)jpeg_turbo_bench(tjtime_fast, mmap.data(), mmap.size(), true);
            }
            if (i % 100 == 0)
                std::cout << "." << std::flush;
        }
        std::cout << std::endl;


        auto imgMPixSize = image.width * image.height * 1e-6;

        std::cout << "image  = " << image.width << "x" << image.height << " (" << std::fixed << std::setprecision(1) << imgMPixSize << " MPix)" << std::endl;
        if (turbo_ok)
        {
            print_stat("jpeg-turbo", tjtime.elapsed(), imgMPixSize, image.size, times);
            print_stat("jpeg-turbo (fast IDCT)", tjtime_fast.elapsed(), imgMPixSize, image.size, times);
            std::cout << "ratio = " << std::fixed << std::setprecision(2) << (tjtime.elapsed() / tjtime_fast.elapsed()) << "x" << std::endl;

        }
        print_stat("nanojpeg", njtime.elapsed(), imgMPixSize, image.size, times);
        if (turbo_ok)
            std::cout << "ratio = " << std::fixed << std::setprecision(2) << (tjtime.elapsed() / njtime.elapsed()) << "x" << std::endl;

        auto out_filename = std::string(argv[1]) + ".bmp";

        std::cout << image.components.size() << " components" << std::endl;
        int comp_nb = std::min<int>(3, image.components.size());
        std::vector<uint8_t> rgb(image.width * image.height * comp_nb);
        convert_time.start();
        auto rgb_out = rgb.data();
        convert(image, [&comp = image.components](int comp_n, int x, int y)
                { return (int)comp[comp_n].pixels[y * comp[comp_n].stride + x]; }, [&rgb_out](int x, int y, auto &&...args)
                { ((*rgb_out++ = args), ...); });
        convert_time.stop();
        std::cout << "yuv to rgb time = " << convert_time.elapsed_str() << std::endl;
        stbi_write_bmp(out_filename.c_str(), image.width, image.height, comp_nb, rgb.data());
    }
    catch (const std::exception &e)
    {
        std::cout << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }
}