# YUV Image Compression with Vector Quantization
A Java program for compressing and decompressing images in the YUV color space using vector quantization.
## Features

- Converts RGB images to YUV color space
- Subsamples U and V channels for enhanced compression
- Generates codebooks for Y, U, and V channels from a training dataset
- Compresses images using vector quantization with 2x2 blocks
- Decompresses images using pre-generated codebooks
- Calculates Mean Squared Error (MSE) and compression ratio
- Supports PNG, JPG, and JPEG input formats; outputs decompressed images as JPG
- Multithreaded processing for codebook generation

## Requirements

1. Java 8 or higher
2. Input images in PNG, JPG, or JPEG format
3. Training images in a specified folder for codebook generation


## Usage

1. Clone the repository:

   ```
    git clone https://github.com/mariamelghandoor/YUV-Image-Compression-with-Vector-Quantization.git


2. Navigate to the project directory:

   ```
    cd yuv-image-compression


3. Compile and run:

    ```
    javac App.java
    java App


4. Enter the image file name (e.g., image.jpg) from the Data\Test\ folder when prompted.

5. Ensure training images are in Data\Train\ and codebooks are in Data\Codebooks\.


## Notes

- The program uses a fixed block size of 2x2 pixels and a codebook size of 256 entries.
- Paths are hardcoded for a specific user directory (C:\Users\mariam\...); update them in the source code as needed.
- Compression ratio is calculated as a percentage of the original size, accounting for YUV subsampling.
- The codebook generation section is commented out in the main method; uncomment to regenerate codebooks.
- Subsampling reduces U and V channels to half the dimensions of the Y channel (4:2:0 chroma subsampling).

## License
MIT License

