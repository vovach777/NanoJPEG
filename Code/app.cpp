#include <iostream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <mio/mmap.hpp>
#include <algorithm>
#include "profiling.hpp"
#include "converter.hpp"
#include "nanojpeg.hpp"

using profiling::StopWatch;

auto nanojpeg_bench(StopWatch &bench, const uint8_t *buf, size_t size)
{
    bench.start();
    auto frame = nanojpeg::decode(buf, size);
    bench.stop();
    bench.elapsed_total -= frame.allocation_time;
    return frame;
}

void nanojpeg_motion_bench(StopWatch &bench, const uint8_t *buf, size_t size, nanojpeg::nj_result &frame)
{
    bench.start();
    nanojpeg::decode(buf, size, frame);
    bench.elapsed_total -= frame.allocation_time;
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
        std::vector<uint8_t> jpeg_vector(mmap.begin(),mmap.end());
        StopWatch njtime{};
        auto file_extenstion = ext(argv[1]);
        bool is_motion = file_extenstion == "mjpeg" ||  file_extenstion == "mjpg";
        if ( is_motion ) {
            std::cout << "Motion JPEG detected" << std::endl;
            StopWatch motion_time{};

            const uint8_t * pos = jpeg_vector.data();
            size_t    size = jpeg_vector.size();
            int times = 0; // number of frames decoded
            std::ofstream fileyuv_file( std::string(argv[1]) + ".y4m", std::ios::binary | std::ios::out | std::ios::trunc);
            nanojpeg::nj_result frame{};

            try
            {
            while (size > 0)
            {
                nanojpeg_motion_bench(motion_time, pos, size,frame);
                pos   += frame.size;
                size  -= frame.size;
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
            }
            catch (const nanojpeg::nj_exception & e) {
                std::cerr << std::endl << "truncated mjpeg stream: " <<  e.what() << " NJ_ERROR: " << e.value <<  std::endl;
            }
            std::cout << std::endl;
            auto imgMPixSize = frame.width * frame.height * 1e-6;
            std::cout << "stream  = " << frame.width << "x" << frame.height << " (" << std::fixed << std::setprecision(1) << imgMPixSize << " MPix)" << std::endl;

            print_stat("nanojpeg motion", motion_time.elapsed(), imgMPixSize,  mmap.size(), times);
            return 0;
        }

        auto image = nanojpeg::decode(jpeg_vector.data(), jpeg_vector.size());//load memmap into nanojpeg image
        if (image.yuv_format == 0) {
            std::cerr << "warining: not standard yuv format!" << std::endl;
        } else {
            std::cout << "* YUV" << image.yuv_format << std::endl;
            std::cout << "* " << image.width << "x" << image.height << std::endl;
            std::cout << "* " << std::fixed << std::setprecision(1) << (image.width*image.height*image.planes.size() / double(image.size)) << ":1" << std::endl;
        }

        int times = std::clamp(1024 * 1024 * 256.0 / image.size / 3,1.,10000.);

        for (int i = 0; i < times; i++)
        {

            (void)nanojpeg_bench(njtime, jpeg_vector.data(), jpeg_vector.size());

            if (i % 100 == 0)
                std::cout << "." << std::flush;
        }
        std::cout << std::endl;


        auto imgMPixSize = image.width * image.height * 1e-6;

        std::cout << "image  = " << image.width << "x" << image.height << " (" << std::fixed << std::setprecision(1) << imgMPixSize << " MPix)" << std::endl;
        print_stat("nanojpeg", njtime.elapsed(), imgMPixSize, image.size * times, times);

        if ( image.planes.size() == 3 && !image.is_ycck)
        {
            auto out_filename = std::string(argv[1]) + ".ppm";
            std::cout << "Creating " << out_filename << std::flush;

            std::ofstream ppm;
            ppm.exceptions(std::ifstream::failbit | std::ifstream::badbit);
            ppm.open(out_filename,  std::ios::binary | std::ios::out);

            std::cout << std::endl << "YUV -> RGB" << std::flush;
            int comp_nb = std::min<int>(3, image.planes.size());
            ppm << "P3\n";
            ppm << image.width << " " << image.height << "\n";
            ppm << "255\n";

            convert<false,3>([&comp = image.planes](int comp_n, int x, int y)
                    {   x >>= comp[comp_n].chroma_w_log2;
                        y >>= comp[comp_n].chroma_h_log2;
                        return (int)comp[comp_n].pixels[ y * comp[comp_n].stride + x]; },

                    [&ppm](int x, int y, int r, int g, int b)
                    {
                        ppm << r << " " << g << " " << b << "\n";
                        if (y % 256==0 && x==0)
                            std::cout << ".";
                    },image.width, image.height);
            std::cout << std::endl;
        }
    }
    catch (const std::exception &e)
    {
        std::cout << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }
}