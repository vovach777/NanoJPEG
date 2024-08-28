#include <iostream>
#include <sstream>
#include "nanojpeg.hpp"
#include "mio.hpp"
#include <algorithm>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "turbojpeg.h"

#ifndef stbi__float2fixed
#define stbi__float2fixed(x)  (((int) ((x) * 4096.0f + 0.5f)) << 8)
#endif
template <typename T>
constexpr inline void YCbCr_to_RGB(int y, int cb, int cr,  T& r_, T& g_, T& b_)
{
      int y_fixed = (y << 20) + (1<<19); // rounding
      int r{},g{},b{};
      cr -= 128;
      cb -= 128;
      r = y_fixed +  cr* stbi__float2fixed(1.40200f);
      g = y_fixed + (cr*-stbi__float2fixed(0.71414f)) + ((cb*-stbi__float2fixed(0.34414f)) & 0xffff0000);
      b = y_fixed                                     +   cb* stbi__float2fixed(1.77200f);
      r >>= 20;
      g >>= 20;
      b >>= 20;
      if ((unsigned) r > 255) { if (r < 0) r = 0; else r = 255; }
      if ((unsigned) g > 255) { if (g < 0) g = 0; else g = 255; }
      if ((unsigned) b > 255) { if (b < 0) b = 0; else b = 255; }
      r_ = r;
      g_ = g;
      b_ = b;
}


template <int planes_nb,typename YUV_, typename RGB_, typename COMP>
inline void convert(YUV_ && yuv, RGB_ && rgb, COMP && comp) {

    int width = comp[0].width;
    int height = comp[0].height;
    for (int h=0; h < height; ++h)
    for (int x=0; x < width;  ++x)
    {
        if constexpr (planes_nb == 1) {
            rgb( x, h,  std::clamp( yuv(0,x,h), 0,255) );
        } else
        if constexpr (planes_nb == 2) {
           rgb(x,h, std::clamp( yuv(0,x,h), 0,255), std::clamp( yuv(1,x,h), 0, 255));
        } else {
            auto y = yuv(0,x,h);
            auto u = yuv(1,x >> comp[1].chroma_w_log2, h >> comp[1].chroma_h_log2);
            auto v = yuv(2,x >> comp[2].chroma_w_log2, h >> comp[2].chroma_h_log2);
            int r,g,b;
            YCbCr_to_RGB(y,u,v, r,g,b);
            if constexpr (planes_nb == 4)
                rgb(x,h,r,g,b,std::clamp( yuv(3,x,h), 0,255));
            else
                rgb(x,h,r,g,b);
        }
    }
}


#include <chrono>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace profiling {
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;

struct StopWatch {
  high_resolution_clock::time_point start_point{};
  high_resolution_clock::time_point end_point{};

  void start() {
    start_point = high_resolution_clock::now();
  }

  auto elipsed() const {

    return  end_point < start_point ? duration<double>(high_resolution_clock::now() - start_point).count() :  duration<double>(end_point - start_point).count();

  }
  auto stop() {
      end_point = high_resolution_clock::now();

  }
  std::string elipsed_str() const {
        return std::string( (std::stringstream() << std::fixed << std::setprecision(9) << elipsed()).str() );
    };
};
}
using profiling::StopWatch;


     struct comp {
        std::vector<uint8_t> data{};
        int chroma_h_log2{}, chroma_w_log2{};
        void setup(int subsamp, int comp) {
            switch (comp) {
                    case 0: chroma_h_log2 = 0; chroma_w_log2 = 0; break;
                    case 1:
                    case 2:
                        switch (subsamp) {
                            case TJSAMP_444: chroma_h_log2 = 0; chroma_w_log2 = 0; break;
                            case TJSAMP_420: chroma_h_log2 = 1; chroma_w_log2 = 1; break;
                            case TJSAMP_422: chroma_h_log2 = 0; chroma_w_log2 = 1; break;
                        }
            }


        }
     };


void jpeg_turbo_bench(const uint8_t * buf, size_t size ) {

    //StopWatch turbo_time{};
    //turbo_time.start();
     auto handle = tjInitDecompress();
     int width, height, jpegSubsamp, jpegColorspace;
     tjDecompressHeader3(handle, buf, size, &width, &height, &jpegSubsamp, &jpegColorspace);
     //tjPlaneSizeYUV(0, width,0, height, jpegSubsamp);
     comp comps[3]{};

     comps[0].data.resize( tjPlaneSizeYUV(0, width,0, height, jpegSubsamp ) );
     comps[0].setup(jpegSubsamp, 0);

     comps[1].data.resize( tjPlaneSizeYUV(1, width,0, height, jpegSubsamp ) );
     comps[1].setup(jpegSubsamp, 1);

     comps[2].data.resize( tjPlaneSizeYUV(2, width,0, height, jpegSubsamp ) );
     comps[2].setup(jpegSubsamp, 2);

     uint8_t * yuv_planes[3] = { comps[0].data.data(), comps[1].data.data(), comps[2].data.data() };
     int succ = tjDecompressToYUVPlanes(handle, buf, size, (uint8_t **) yuv_planes, 0, nullptr, 0, 0);

     if (succ != 0)
     {
        auto err = std::string(tjGetErrorStr2(handle));
        tjDestroy(handle);
        throw std::runtime_error(err);
     }
     tjDestroy(handle);
     //turbo_time.stop();
     //return turbo_time.elipsed();
}

void nanojpeg_bench(const uint8_t * buf, size_t size ) {
    NanoJpeg decoder;
    decoder.njDecode(buf, size);
}



inline void print_stat(std::string title, double elipsed, double imgMPixSize, size_t size, int times)
{
    std::cout << "** " << title << "  **" << std::endl;
    std::cout << "time   = " << std::fixed << std::setprecision(9) << (elipsed/ times) << " seconds" << std::endl;
    std::cout << "images = " << std::fixed << std::setprecision(3) << (times / elipsed) << " fps" << std::endl;
    std::cout << "bytes  = " << std::fixed << std::setprecision(2) <<  (times * size / ( 1024.0 * 1024.0 ) / elipsed ) << " MB/s" << std::endl;
    std::cout << "pixels = " << std::fixed << std::setprecision(2) <<  (imgMPixSize * double( times ) / elipsed) << " MPix/s" << std::endl;
    std::cout << std::endl;
}

int main(int argc, char ** argv) {

 try {
    std::error_code error;
    mio::ummap_source mmap = mio::make_mmap<mio::ummap_source>(argv[1], 0, 0, error);
    if (error)
    {
        throw std::runtime_error(error.message());
    }
    StopWatch decode_time, convert_time;

    int times = std::min( std::ceil( 1024*1024*512.0 / mmap.size() ), 100.);
    std::cout << "bencmarking turbo" << std::endl;
    decode_time.start();
    for (int i = 0; i < times; i++)
    {
        jpeg_turbo_bench(mmap.data(),mmap.size());
    }
    decode_time.stop();
    std::cout << "bencmarking nano" << std::endl;

    auto turbo_elipsed = decode_time.elipsed();


    decode_time.start();
    for (int i = 0; i < times; i++)
    {
        nanojpeg_bench(mmap.data(),mmap.size());
    }
    decode_time.stop();


    auto elipsed = decode_time.elipsed();

    NanoJpeg decoder;
    decoder.njDecode(mmap.data(),mmap.size());

    auto imgMPixSize = decoder.njGetWidth() * decoder.njGetHeight() *  1e-6;


    std::cout << "image  = " << decoder.njGetWidth() << "x" << decoder.njGetHeight() << " (" << std::fixed << std::setprecision(3) << imgMPixSize << " MPix)" << std::endl;
    print_stat("jpeg-turbo",turbo_elipsed,imgMPixSize,mmap.size(), times);
    print_stat("nanojpeg",elipsed,imgMPixSize,mmap.size(), times);
    std::cout << "ratio = " << std::fixed << std::setprecision(2) << (turbo_elipsed / elipsed) << "x" << std::endl;

    auto out_filename = std::string( argv[1] ) + ".bmp";
    auto& comp = decoder.njGetComponents();
    if (comp.size() == 3) {
        std::cout << "3 components" << std::endl;
        std::vector<uint8_t> rgb( comp[0].width * comp[0].height * 3);
        convert_time.start();
         convert<3>([&comp](int comp_n, int x, int y ) {
            return (int)comp[comp_n].pixels[y*comp[comp_n].stride+x];
         },[width=comp[0].width,&rgb](int x, int y, auto r, auto g, auto b){
            auto pos =   rgb.begin()  + (y * width + x)*3;
            *pos++ = r;
            *pos++ = g;
            *pos++ = b;
         }, comp);
         convert_time.stop();
         std::cout << "yuv to rgb time = " << convert_time.elipsed_str() << std::endl;
        stbi_write_bmp( out_filename.c_str(),decoder.njGetWidth(), decoder.njGetHeight(), 3, rgb.data() );
    }
 } catch(const std::exception & e) {
    std::cerr << e.what() << std::endl;
    return 1;
 }

}