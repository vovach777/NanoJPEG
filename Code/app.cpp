#include <iostream>
#include <sstream>
#include <tuple>
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



template <bool is_ycck, int planes_nb,typename YUV_, typename RGB_, typename COMP>
inline void convert(YUV_ && yuv, RGB_ && rgb, COMP && comp, int width, int height) {

    for (int h=0; h < height; ++h)
    for (int x=0; x < width;  ++x)
    {
        if constexpr (planes_nb == 1) {
            rgb( x, h,  std::clamp( yuv(0,x,h), 0,255) );
        } else
        if constexpr (planes_nb == 2) {
           rgb(x,h, std::clamp( yuv(0,x,h), 0,255), std::clamp( yuv(1,x,h), 0, 255));
        } else
        if constexpr (planes_nb == 3) {
            auto y = yuv(0,x,h);
            auto u = yuv(1,x >> comp[1].chroma_w_log2, h >> comp[1].chroma_h_log2);
            auto v = yuv(2,x >> comp[2].chroma_w_log2, h >> comp[2].chroma_h_log2);
            int r,g,b;
            YCbCr_to_RGB(y,u,v, r,g,b);
            rgb(x,h,r,g,b);
        } else {
            const auto c =  yuv(0,x >> comp[0].chroma_w_log2, h >> comp[0].chroma_h_log2);
            const auto m = yuv(1,x >> comp[1].chroma_w_log2, h >> comp[1].chroma_h_log2);
            const auto y = yuv(2,x >> comp[2].chroma_w_log2, h >> comp[2].chroma_h_log2);
            const auto k =  yuv(3,x >> comp[3].chroma_w_log2, h >> comp[3].chroma_h_log2);
            if constexpr (is_ycck) {
                int ir,ig,ib;
                YCbCr_to_RGB(c,m,y, ir,ig,ib);
                const auto r = nanojpeg::njClip( (255-ir)*k/255 );
                const auto g = nanojpeg::njClip( (255-ig)*k/255 );
                const auto b = nanojpeg::njClip( (255-ib)*k/255 );
                rgb(x,h,r,g,b);
            } else {
                const auto r = nanojpeg::njClip( c*k/255 );
                const auto g = nanojpeg::njClip( m*k/255 );
                const auto b = nanojpeg::njClip( y*k/255 );
                rgb(x,h,r,g,b);
            }
        }
    }
}

template <typename YUV_, typename RGB_>
inline void convert(nanojpeg::nj_result &image,  YUV_ && yuv, RGB_ && rgb)
{
    switch ( image.components.size() )
    {
        case 1: convert<false,1>( std::forward<YUV_>(yuv), std::forward<RGB_>(rgb), image.components,image.width,image.height); break;
        case 2: convert<false,2>( std::forward<YUV_>(yuv), std::forward<RGB_>(rgb), image.components,image.width,image.height); break;
        case 3: convert<false,3>( std::forward<YUV_>(yuv), std::forward<RGB_>(rgb), image.components,image.width,image.height); break;
        case 4: if (image.is_ycck)
                    convert<true,4>( std::forward<YUV_>(yuv), std::forward<RGB_>(rgb), image.components,image.width,image.height);
                else
                    convert<false,4>( std::forward<YUV_>(yuv), std::forward<RGB_>(rgb), image.components,image.width,image.height);
                break;
        default:
            throw std::runtime_error("invalid number of planes");
    }


}


#include <chrono>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <optional>
namespace profiling {
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;

struct StopWatch {
  std::optional<high_resolution_clock::time_point> start_point{};
  double elipsed_total{};

  inline void startnew() {
    elipsed_total = 0.0;
    start();
  }

  inline void start() {
   if ( !start_point )
        start_point = high_resolution_clock::now();
  }

  inline bool is_running() const { return start_point.has_value(); }

  inline double elipsed() const {
    return elipsed_total + ( (start_point) ?  duration<double>(high_resolution_clock::now() - *start_point).count() : 0 );

  }
  void stop() {
    if ( start_point ) {
        elipsed_total += duration<double>(high_resolution_clock::now() - *start_point).count();
        start_point.reset();
    }

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

    (void)nanojpeg::decode(buf, size);
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
    StopWatch njtime{}, tjtime{}, convert_time{};

    int times = std::max( std::ceil( 1024*1024*100.0 / mmap.size() ), 60.);
    std::cout << "benchmarking" << std::flush;

    bool turbo_ok = true;
    try {
    jpeg_turbo_bench(mmap.data(),mmap.size());
    } catch(std::exception & e) {
        std::cout << "turbojpeg failed: " << e.what() << std::endl;
        turbo_ok = false;
    }

    for (int i = 0; i < times; i++)
    {
        njtime.start();
        nanojpeg_bench(mmap.data(),mmap.size());
        njtime.stop();
        if ( turbo_ok ) {
            tjtime.start();
            jpeg_turbo_bench(mmap.data(),mmap.size());
            tjtime.stop();
        }
        if (i % 50 == 0) std::cout << "." << std::flush;
    }
    std::cout << std::endl;





    auto image = nanojpeg::decode(mmap.data(),mmap.size());

    auto imgMPixSize = image.width * image.height *  1e-6;


    std::cout << "image  = " << image.width << "x" << image.height << " (" << std::fixed << std::setprecision(3) << imgMPixSize << " MPix)" << std::endl;
    if (turbo_ok)
        print_stat("jpeg-turbo",tjtime.elipsed(),imgMPixSize,mmap.size(), times);
    print_stat("nanojpeg",njtime.elipsed(),imgMPixSize,mmap.size(), times);
    if (turbo_ok)
        std::cout << "ratio = " << std::fixed << std::setprecision(3) << (tjtime.elipsed() / njtime.elipsed()) << "x" << std::endl;

    auto out_filename = std::string( argv[1] ) + ".bmp";



        std::cout << image.components.size() << " components" << std::endl;
        int comp_nb = std::min<int>(3, image.components.size());
        std::vector<uint8_t> rgb( image.width * image.height * comp_nb);
        convert_time.start();
        auto rgb_out = rgb.data();
         convert(image, [&comp=image.components](int comp_n, int x, int y ) {
            return (int)comp[comp_n].pixels[y*comp[comp_n].stride+x];
         },[&rgb_out](int x, int y, auto&&... args){
           
            ( (*rgb_out++ = args), ...);
         });
         convert_time.stop();
         std::cout << "yuv to rgb time = " << convert_time.elipsed_str() << std::endl;
        stbi_write_bmp( out_filename.c_str(),image.width, image.height, comp_nb, rgb.data() );

 } catch(const std::exception & e) {
    std::cerr << e.what() << std::endl;
    return 1;
 }

}