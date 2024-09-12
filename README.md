# NanoJPEG++

**NanoJPEG++** is a project created with the goal of optimizing the original **NanoJPEG** library to achieve performance close to that of **libjpeg-turbo**. This version of the library is developed in C++ and focuses on delivering high-speed JPEG image decoding. **NanoJPEG++** is a header-only project, making integration into your project straightforward.

## Description

**NanoJPEG++** is a modern, high-performance library for JPEG image decoding. This version includes several enhancements to achieve a significant performance boost compared to the original **NanoJPEG** library:

- **Bitstream Processing**: Reading data from the stream in 32-bit chunks allows simultaneous extraction of Huffman codes and their symbols. This speeds up decoding by reducing the number of shift operations.
- **Optimized Lookup Table**: The Huffman table has been reduced and optimized to use less memory. It is split into DC and AC components with different bit depths. If a direct lookup fails, the number of leading ones is computed, and a quick transition is made to the nearest Huffman code to optimize the search.
- **IDCT Vectorization**: Modern compilers have achieved automatic vectorization of the IDCT (Discrete Cosine Transform) process, eliminating the need for complex vector operation code or assembly inserts.

**NanoJPEG++** is a **header-only** project, which means you do not need to build separate binaries or libraries. Simply include the header files in your project to take advantage of the library's features.

## Changes from NanoJPEG

In the original **NanoJPEG** library, RGB conversion was a significant performance bottleneck. In **NanoJPEG++**, I chose to omit the RGB conversion and focus solely on decoding JPEG components. For color space conversion, you will need to use other specialized libraries.

## Example Usage

Here's how you can use **NanoJPEG++** to decode a JPEG image:

```cpp
#include "nanojpeg.h"

int main() {
    const uint8_t* buf; // JPEG data buffer
    size_t size;        // Buffer size

    // Assume buf and size are already initialized
    auto result = nanojpeg::decode(buf, size);

    // Access the decoding results
    int width = result.width;
    int height = result.height;
    auto& components = result.components;
    bool is_ycck = result.is_ycck;
}
```

## Performance

**NanoJPEG++** shows an average performance improvement of 5% compared to `libjpeg-turbo`, especially on larger files. This makes it an excellent choice for applications requiring high-speed image processing.

- The `app.cpp` file in the project contains usage examples and can be used for performance testing.

## License

**NanoJPEG++** is licensed under the MIT License. See [LICENSE](LICENSE) for more details.

## Contributing

We welcome contributions to the project. If you would like to make improvements, please create a pull request or report issues in the GitHub Issues section.

## Contact

If you have any questions or suggestions, you can reach us at [vovach777@gmail.com](mailto:vovach777@gmail.com).
