#include <iostream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <mio/mmap.hpp>
#include <algorithm>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "turbojpeg.h"
#include "profiling.hpp"
#include "converter.hpp"
#include "nanojpeg.hpp"

using profiling::StopWatch;

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

static auto nanojpeg_bench(StopWatch &bench, const uint8_t *buf, size_t size)
{
    bench.start();
    auto frame = nanojpeg::decode(buf, size);
    bench.stop();
    bench.elapsed_total -= frame.allocation_time;
    return frame;
}

static void nanojpeg_motion_bench(StopWatch &bench, const uint8_t *buf, size_t size, nanojpeg::nj_result &frame)
{
    bench.start();
    nanojpeg::decode(buf, size, frame);
    bench.stop();
}

inline void print_stat(std::string title, double elapsed, double imgMPixSize, size_t streamSize, int times)
{
    std::cout << "** " << title << " **" << std::endl;
    std::cout << "time   = " << std::fixed << std::setprecision(9) << (elapsed / times) << " seconds" << std::endl;
    std::cout << "images = " << std::fixed << std::setprecision(3) << (times / elapsed) << " fps" << std::endl;
    std::cout << "bytes  = " << std::fixed << std::setprecision(2) << (streamSize / (1024.0 * 1024.0) / elapsed) << " MB/s" << std::endl;
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
                pos   += frame.size;
                size  -= frame.size;
                motion_time.stop();
                times+=1;
                if (times == 1) {
                    //header
                    if ( frame.yuv_format != 420 )
                        fileyuv_file << "YUV4MPEG2 W" << frame.width << " H" << frame.height << " F25:1 It A1:1 C" << frame.yuv_format <<  " XYSCSS=" << frame.yuv_format << " XCOLORRANGE=FULL\n";
                    else
                        fileyuv_file << "YUV4MPEG2 W" << frame.width << " H" << frame.height << " F25:1 It A1:1 C420jpeg XYSCSS=420JPEG XCOLORRANGE=FULL\n";

                }
                if (times % 100 == 0)
                    std::cout << "." << std::flush;

                fileyuv_file << "FRAME\n";

                for (const auto& c : frame.planes)
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

        int times = std::clamp(1024 * 1024 * 256.0 / image.size / 3,1.,10000.);
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
            print_stat("jpeg-turbo", tjtime.elapsed(), imgMPixSize, image.size * times, times);
            print_stat("jpeg-turbo (fast IDCT)", tjtime_fast.elapsed(), imgMPixSize, image.size * times, times);
            std::cout << "ratio = " << std::fixed << std::setprecision(2) << (tjtime.elapsed() / tjtime_fast.elapsed()) << "x" << std::endl;

        }
        print_stat("nanojpeg", njtime.elapsed(), imgMPixSize, image.size * times, times);
        if (turbo_ok)
            std::cout << "ratio = " << std::fixed << std::setprecision(2) << (tjtime.elapsed() / njtime.elapsed()) << "x" << std::endl;

        auto out_filename = std::string(argv[1]) + ".bmp";

        std::cout << std::endl << "YUV -> RGB..." << std::flush;
        int comp_nb = std::min<int>(3, image.planes.size());
        std::vector<uint8_t> rgb(image.width * image.height * comp_nb);
        convert_time.start();
        auto rgb_out = rgb.data();
        convert(image.width, image.height, image.planes.size(), image.is_ycck,  [&comp = image.planes](int comp_n, int x, int y)
                {   x >>= comp[comp_n].chroma_w_log2;
                    y >>= comp[comp_n].chroma_h_log2;
                    return (int)comp[comp_n].pixels[ y * comp[comp_n].stride + x]; }, [&rgb_out](int x, int y, auto &&...args)
                { ((*rgb_out++ = args), ...); });
        convert_time.stop();
        std::cout << std::endl << "time = " << convert_time.elapsed_str() << std::endl;
        std::cout << "Creating " << out_filename << std::flush;
        stbi_write_bmp(out_filename.c_str(), image.width, image.height, comp_nb, rgb.data());
        std::cout << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cout << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }
}