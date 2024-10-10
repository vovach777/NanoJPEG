#pragma once
#include <utility>
#include <cstdint>
#include <algorithm>
#include <tuple>
constexpr inline uint8_t converter_u8(float v) {
    return std::clamp<int>( v+0.5f, 0,255);
}
constexpr inline auto YCbCr_to_RGB(int y, int cb, int cr)
{
    const float cr_ = cr - 128;
    const float cb_ = cb - 128;
    return std::make_tuple(
            converter_u8(y + cr_ * 1.40200f),
            converter_u8(y + cr_ * -0.71414f + cb_ * -0.34414f),
            converter_u8(y + cb_ * 1.77200f) );
}

template <bool is_ycck, int planes_nb, typename YUV_, typename RGB_>
inline void convert(YUV_ &&yuv, RGB_ &&rgb, int width, int height)
{

    for (int h = 0; h < height; ++h)
        for (int x = 0; x < width; ++x)
        {
            if constexpr (planes_nb == 1)
            {
                rgb(x, h, yuv(0, x, h));
            }
            else if constexpr (planes_nb == 2)
            {
                rgb(x, h, yuv(0, x, h), yuv(1, x, h));
            }
            else if constexpr (planes_nb == 3)
            {
                auto y = yuv(0, x, h);
                auto u = yuv(1, x, h);
                auto v = yuv(2, x, h);
                auto [r, g, b] = YCbCr_to_RGB(y, u, v);
                rgb(x, h, r, g, b);
            }
            else
            {
                const auto c = yuv(0, x, h);
                const auto m = yuv(1, x, h);
                const auto y = yuv(2, x, h);
                const auto k = yuv(3, x, h);
                if constexpr (is_ycck)
                {
                    auto [ir, ig, ib] = YCbCr_to_RGB(c, m, y);
                    const auto r = converter_u8((255 - ir) * k / 255);
                    const auto g = converter_u8((255 - ig) * k / 255);
                    const auto b = converter_u8((255 - ib) * k / 255);
                    rgb(x, h, r, g, b);
                }
                else
                {
                    const auto r = converter_u8(c * k / 255);
                    const auto g = converter_u8(m * k / 255);
                    const auto b = converter_u8(y * k / 255);
                    rgb(x, h, r, g, b);
                }
            }
        }
}

template <typename YUV_, typename RGB_>
inline void convert(int width, int height,int planes, bool is_ycck, YUV_ &&yuv, RGB_ &&rgb)
{
    switch (planes)
    {
    case 1:
        convert<false, 1>(std::forward<YUV_>(yuv), std::forward<RGB_>(rgb), width, height);
        break;
    case 2:
        convert<false, 2>(std::forward<YUV_>(yuv), std::forward<RGB_>(rgb), width, height);
        break;
    case 3:
        convert<false, 3>(std::forward<YUV_>(yuv), std::forward<RGB_>(rgb), width, height);
        break;
    case 4:
        if (is_ycck)
            convert<true, 4>(std::forward<YUV_>(yuv), std::forward<RGB_>(rgb), width, height);
        else
            convert<false, 4>(std::forward<YUV_>(yuv), std::forward<RGB_>(rgb), width, height);
        break;
    default:
        throw std::runtime_error("invalid number of planes");
    }
}
