import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.Scanner;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;

public class App {

    public static List<List<int[]>> collectAllYUVBlocks(String folderPath, int blockSize) throws IOException {
        List<int[]> allYBlocks = new ArrayList<>();
        List<int[]> allUBlocks = new ArrayList<>();
        List<int[]> allVBlocks = new ArrayList<>();
        
        File folder = new File(folderPath);
        if (!folder.exists() || !folder.isDirectory()) {
            throw new IOException("Invalid directory path: " + folderPath);
        }
        
        File[] files = folder.listFiles((dir, name) -> 
            name.toLowerCase().endsWith(".jpg") || 
            name.toLowerCase().endsWith(".png") || 
            name.toLowerCase().endsWith(".jpeg")
        );
        
        if (files == null || files.length == 0) {
            throw new IOException("No image files found in directory: " + folderPath);
        }
        
        File[] limitedFiles = Arrays.copyOf(files, files.length);
        
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);
        List<List<int[]>>[] results = new List[files.length];
        
        for (int i = 0; i < files.length; i++) {
            final File file = limitedFiles[i];
            final int index = i;
            executor.submit(() -> {
                try {
                    BufferedImage image = ImageIO.read(file);
                    List<List<int[]>> yuvBlocks = getYUVBlocks(image, blockSize);
                    results[index] = yuvBlocks;
                } catch (IOException e) {
                    System.err.println("Error processing file " + file.getName() + ": " + e.getMessage());
                }
            });
        }
        
        executor.shutdown();
        try {
            executor.awaitTermination(1, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            throw new IOException("Image processing interrupted", e);
        }
        
        for (List<List<int[]>> yuvBlocks : results) {
            if (yuvBlocks != null) {
                allYBlocks.addAll(yuvBlocks.get(0));  
                allUBlocks.addAll(yuvBlocks.get(1));
                allVBlocks.addAll(yuvBlocks.get(2));  
            }
        }
        
        List<List<int[]>> result = new ArrayList<>();
        result.add(allYBlocks);
        result.add(allUBlocks);
        result.add(allVBlocks);
        return result;
    }
    
    public static List<List<int[]>> getYUVBlocks(BufferedImage image, int blockSize) throws IOException {
        List<int[]> yBlocks = new ArrayList<>();
        List<int[]> uBlocks = new ArrayList<>();
        List<int[]> vBlocks = new ArrayList<>();
        
        int width = image.getWidth();
        int height = image.getHeight();
        
        double[][] yMatrix = new double[height][width];
        double[][] uMatrix = new double[height][width];
        double[][] vMatrix = new double[height][width];
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int rgb = image.getRGB(x, y);
                int r = (rgb >> 16) & 0xFF;
                int g = (rgb >> 8) & 0xFF;
                int b = rgb & 0xFF;
                
                int yVal = (int) Math.round(0.299 * r + 0.587 * g + 0.114 * b);
                int uVal = (int) Math.round(-0.147 * r - 0.289 * g + 0.436 * b);
                int vVal = (int) Math.round(0.615 * r - 0.515 * g - 0.100 * b);
                
                yMatrix[y][x] = Math.max(0, Math.min(255, yVal));
                uMatrix[y][x] = uVal + 128; 
                vMatrix[y][x] = vVal + 128;
            }
        }
        
        double[][] subsampledU = subsampleChannel(uMatrix);
        double[][] subsampledV = subsampleChannel(vMatrix);
        
        for (int y = 0; y <= height - blockSize; y += blockSize) {
            for (int x = 0; x <= width - blockSize; x += blockSize) {
                int[] yBlock = new int[blockSize * blockSize];
                int index = 0;
                for (int dy = 0; dy < blockSize; dy++) {
                    for (int dx = 0; dx < blockSize; dx++) {
                        yBlock[index++] = (int) yMatrix[y + dy][x + dx];
                    }
                }
                yBlocks.add(yBlock);
            }
        }
        
        for (int y = 0; y <= subsampledU.length - blockSize; y += blockSize) {
            for (int x = 0; x <= subsampledU[0].length - blockSize; x += blockSize) {
                int[] uBlock = new int[blockSize * blockSize];
                int index = 0;
                for (int dy = 0; dy < blockSize; dy++) {
                    for (int dx = 0; dx < blockSize; dx++) {
                        uBlock[index++] = (int) subsampledU[y + dy][x + dx];
                    }
                }
                uBlocks.add(uBlock);
            }
        }
        
        for (int y = 0; y <= subsampledV.length - blockSize; y += blockSize) {
            for (int x = 0; x <= subsampledV[0].length - blockSize; x += blockSize) {
                int[] vBlock = new int[blockSize * blockSize];
                int index = 0;
                for (int dy = 0; dy < blockSize; dy++) {
                    for (int dx = 0; dx < blockSize; dx++) {
                        vBlock[index++] = (int) subsampledV[y + dy][x + dx];
                    }
                }
                vBlocks.add(vBlock);
            }
        }
        
        List<List<int[]>> result = new ArrayList<>();
        result.add(yBlocks);
        result.add(uBlocks);
        result.add(vBlocks);
        
        return result;
    }

    public static double[][] subsampleChannel(double[][] channel) {
        if (channel == null || channel.length == 0 || channel[0].length == 0) {
            throw new IllegalArgumentException("Invalid channel matrix");
        }

        int rows = channel.length;
        int cols = channel[0].length;

        if (rows % 2 != 0 || cols % 2 != 0) {
            throw new IllegalArgumentException("Channel dimensions must be even for subsampling");
        }

        int newRows = rows / 2;
        int newCols = cols / 2;
        double[][] subsampled = new double[newRows][newCols];

        for (int i = 0; i < newRows; i++) {
            for (int j = 0; j < newCols; j++) {
                subsampled[i][j] = channel[i * 2][j * 2];
            }
        }

        return subsampled;
    }

    public static double[][] upsampleChannel(double[][] channel, int targetRows, int targetCols) {
        if (channel == null || channel.length == 0 || channel[0].length == 0) {
            throw new IllegalArgumentException("Invalid channel matrix");
        }

        int rows = channel.length;
        int cols = channel[0].length;

        if (targetRows < rows || targetCols < cols) {
            throw new IllegalArgumentException("Target dimensions must be larger than source dimensions");
        }

        double[][] upsampled = new double[targetRows][targetCols];

        for (int i = 0; i < targetRows; i++) {
            for (int j = 0; j < targetCols; j++) {
                int srcI = Math.min(i / 2, rows - 1);
                int srcJ = Math.min(j / 2, cols - 1);
                upsampled[i][j] = channel[srcI][srcJ];
            }
        }

        return upsampled;
    }

    public static double[] computeAverageVector(List<int[]> blocks) {
        if (blocks == null || blocks.isEmpty()) {
            return new double[0];
        }
        int vectorLength = blocks.get(0).length;
        double[] avgVector = new double[vectorLength];
        for (int[] block : blocks) {
            for (int i = 0; i < vectorLength; i++) {
                avgVector[i] += block[i];
            }
        }
        int blockCount = blocks.size();
        for (int i = 0; i < vectorLength; i++) {
            avgVector[i] /= blockCount;
        }
        return avgVector;
    }

    public static double Compute_Distance(int[] vector1, double[] vector2) {
        double distance = 0;
        for (int i = 0; i < vector1.length; i++) {
            double diff = vector1[i] - vector2[i];
            distance += diff * diff;
        }
        return Math.sqrt(distance);
    }

    public static double[][] splitVector(double[] vector) {
        int length = vector.length;
        double[] plusHalf = new double[length];
        double[] minusHalf = new double[length];
        for (int i = 0; i < length; i++) {
            plusHalf[i] = vector[i] + 5.0;
            minusHalf[i] = vector[i] - 5.0;
        }
        return new double[][]{plusHalf, minusHalf};
    }

    public static HashMap<String, double[]> Create_Code_Book(List<int[]> blocks, int codeBookSize) {
        if (blocks == null || blocks.isEmpty()) {
            return new HashMap<>();
        }
        List<double[]> startvector = new ArrayList<>();
        startvector.add(computeAverageVector(blocks));
        while (startvector.size() < codeBookSize) {
            List<double[]> newvector = new ArrayList<>();
            for (double[] vector : startvector) {
                double[][] split = splitVector(vector);
                double[] plusHalf = split[0];
                double[] minusHalf = split[1];
                List<int[]> plusHalfBlocks = new ArrayList<>();
                List<int[]> minusHalfBlocks = new ArrayList<>();
                for (int[] block : blocks) {
                    double distToPlus = Compute_Distance(block, plusHalf);
                    double distToMinus = Compute_Distance(block, minusHalf);
                    if (distToPlus < distToMinus) {
                        plusHalfBlocks.add(block);
                    } else {
                        minusHalfBlocks.add(block);
                    }
                }
                double[] plusAvg = plusHalfBlocks.isEmpty() ? plusHalf : computeAverageVector(plusHalfBlocks);
                double[] minusAvg = minusHalfBlocks.isEmpty() ? minusHalf : computeAverageVector(minusHalfBlocks);
                newvector.add(plusAvg);
                newvector.add(minusAvg);
            }
            startvector = newvector;
        }
        boolean changed;
        HashMap<Integer, List<int[]>> assignments = new HashMap<>();
        Random random = new Random();
        int iterations = 0;
        do {
            changed = false;
            assignments.clear();
            for (int i = 0; i < startvector.size(); i++) {
                assignments.put(i, new ArrayList<>());
            }
            for (int[] block : blocks) {
                int bestIndex = 0;
                double bestDist = Double.MAX_VALUE;
                for (int i = 0; i < startvector.size(); i++) {
                    double dist = Compute_Distance(block, startvector.get(i));
                    if (dist < bestDist) {
                        bestDist = dist;
                        bestIndex = i;
                    }
                }
                assignments.get(bestIndex).add(block);
            }
            for (int i = 0; i < startvector.size(); i++) {
                if (assignments.get(i).isEmpty() && !blocks.isEmpty()) {
                    assignments.get(i).add(blocks.get(random.nextInt(blocks.size())));
                    changed = true;
                }
            }
            for (int i = 0; i < startvector.size(); i++) {
                List<int[]> group = assignments.get(i);
                if (!group.isEmpty()) {
                    double[] newvector = computeAverageVector(group);
                    if (!Arrays.equals(startvector.get(i), newvector)) {
                        startvector.set(i, newvector);
                        changed = true;
                    }
                }
            }
            iterations++;
        } while (changed && iterations < 100);
        HashMap<String, double[]> codeBook = new HashMap<>();
        int keyLength = (int) Math.ceil(Math.log(codeBookSize) / Math.log(2));
        for (int i = 0; i < startvector.size(); i++) {
            String binaryCode = String.format("%" + keyLength + "s", Integer.toBinaryString(i)).replace(' ', '0');
            codeBook.put(binaryCode, startvector.get(i));
        }
        return codeBook;
    }

    public static void saveCodeBook(HashMap<String, double[]> codeBook, String filePath) throws IOException {
        try (ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(filePath))) {
            oos.writeObject(codeBook);
        }
    }

    public static HashMap<String, double[]> loadCodeBook(String filePath) throws IOException, ClassNotFoundException {
        try (ObjectInputStream ois = new ObjectInputStream(new FileInputStream(filePath))) {
            return (HashMap<String, double[]>) ois.readObject();
        }
    }

    public static void saveCompressedBlocks(List<String> compressedBlocks, String filePath) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            for (String code : compressedBlocks) {
                writer.write(code);
                writer.newLine();
            }
        }
    }

    public static List<String> loadCompressedBlocks(String filePath) throws IOException {
        List<String> compressedBlocks = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = reader.readLine()) != null) {
                compressedBlocks.add(line);
            }
        }
        return compressedBlocks;
    }

    public static List<List<String>> Compress_test_img(BufferedImage image, int blockSize) throws IOException, ClassNotFoundException {
        String codebookDir = "C:\\Users\\mariam\\OneDrive\\Desktop\\vs code\\Java\\final_project\\Data\\Codebooks";
        String yCodebookPath = codebookDir + "\\y_codebook.ser";
        String uCodebookPath = codebookDir + "\\u_codebook.ser";
        String vCodebookPath = codebookDir + "\\v_codebook.ser";
        
        List<List<int[]>> yuvBlocks = getYUVBlocks(image, blockSize);
        List<int[]> yBlocks = yuvBlocks.get(0);
        List<int[]> uBlocks = yuvBlocks.get(1);
        List<int[]> vBlocks = yuvBlocks.get(2);
        
        HashMap<String, double[]> yCodeBook = loadCodeBook(yCodebookPath);
        HashMap<String, double[]> uCodeBook = loadCodeBook(uCodebookPath);
        HashMap<String, double[]> vCodeBook = loadCodeBook(vCodebookPath);
        
        List<String> compressedYBlocks = new ArrayList<>();
        List<String> compressedUBlocks = new ArrayList<>();
        List<String> compressedVBlocks = new ArrayList<>();
        
        for (int[] block : yBlocks) {
            String bestCode = "";
            double minDist = Double.MAX_VALUE;
            for (String code : yCodeBook.keySet()) {
                double dist = Compute_Distance(block, yCodeBook.get(code));
                if (dist < minDist) {
                    minDist = dist;
                    bestCode = code;
                }
            }
            compressedYBlocks.add(bestCode);
        }
        
        for (int[] block : uBlocks) {
            String bestCode = "";
            double minDist = Double.MAX_VALUE;
            for (String code : uCodeBook.keySet()) {
                double dist = Compute_Distance(block, uCodeBook.get(code));
                if (dist < minDist) {
                    minDist = dist;
                    bestCode = code;
                }
            }
            compressedUBlocks.add(bestCode);
        }
        
        for (int[] block : vBlocks) {
            String bestCode = "";
            double minDist = Double.MAX_VALUE;
            for (String code : vCodeBook.keySet()) {
                double dist = Compute_Distance(block, vCodeBook.get(code));
                if (dist < minDist) {
                    minDist = dist;
                    bestCode = code;
                }
            }
            compressedVBlocks.add(bestCode);
        }
        
        List<List<String>> compressedChannels = new ArrayList<>();
        compressedChannels.add(compressedYBlocks);
        compressedChannels.add(compressedUBlocks);
        compressedChannels.add(compressedVBlocks);
        
        return compressedChannels;
    }

    public static BufferedImage Decompress_test_img(String compressedYPath, String compressedUPath, String compressedVPath, int blockSize, int width, int height) throws IOException, ClassNotFoundException {
        String codebookDir = "C:\\Users\\mariam\\OneDrive\\Desktop\\vs code\\Java\\final_project\\Data\\Codebooks";
        String yCodebookPath = codebookDir + "\\y_codebook.ser";
        String uCodebookPath = codebookDir + "\\u_codebook.ser";
        String vCodebookPath = codebookDir + "\\v_codebook.ser";
        
        HashMap<String, double[]> yCodeBook = loadCodeBook(yCodebookPath);
        HashMap<String, double[]> uCodeBook = loadCodeBook(uCodebookPath);
        HashMap<String, double[]> vCodeBook = loadCodeBook(vCodebookPath);
        
        List<String> compressedYBlocks = loadCompressedBlocks(compressedYPath);
        List<String> compressedUBlocks = loadCompressedBlocks(compressedUPath);
        List<String> compressedVBlocks = loadCompressedBlocks(compressedVPath);
        
        BufferedImage decompressedImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        
        // int blocksX = (int) Math.ceil((double) width / blockSize);
        // int blocksY = (int) Math.ceil((double) height / blockSize);
        
        double[][] yMatrix = new double[height][width];
        double[][] uMatrix = new double[height / 2][width / 2];
        double[][] vMatrix = new double[height / 2][width / 2];
        
        // Decompress Y 
        int blockIndex = 0;
        for (int y = 0; y <= height - blockSize && blockIndex < compressedYBlocks.size(); y += blockSize) {
            for (int x = 0; x <= width - blockSize && blockIndex < compressedYBlocks.size(); x += blockSize) {
                double[] yVector = yCodeBook.get(compressedYBlocks.get(blockIndex));
                if (yVector == null) {
                    throw new IOException("Invalid Y code at index " + blockIndex);
                }
                int vectorIndex = 0;
                for (int dy = 0; dy < blockSize && y + dy < height; dy++) {
                    for (int dx = 0; dx < blockSize && x + dx < width; dx++) {
                        yMatrix[y + dy][x + dx] = Math.min(255, Math.max(0, yVector[vectorIndex++]));
                    }
                }
                blockIndex++;
            }
        }
        
        // Decompress U 
        blockIndex = 0;
        for (int y = 0; y <= (height / 2) - blockSize && blockIndex < compressedUBlocks.size(); y += blockSize) {
            for (int x = 0; x <= (width / 2) - blockSize && blockIndex < compressedUBlocks.size(); x += blockSize) {
                double[] uVector = uCodeBook.get(compressedUBlocks.get(blockIndex));
                if (uVector == null) {
                    throw new IOException("Invalid U code at index " + blockIndex);
                }
                int vectorIndex = 0;
                for (int dy = 0; dy < blockSize && y + dy < height / 2; dy++) {
                    for (int dx = 0; dx < blockSize && x + dx < width / 2; dx++) {
                        uMatrix[y + dy][x + dx] = Math.min(255, Math.max(0, uVector[vectorIndex++]));
                    }
                }
                blockIndex++;
            }
        }
        
        // Decompress V 
        blockIndex = 0;
        for (int y = 0; y <= (height / 2) - blockSize && blockIndex < compressedVBlocks.size(); y += blockSize) {
            for (int x = 0; x <= (width / 2) - blockSize && blockIndex < compressedVBlocks.size(); x += blockSize) {
                double[] vVector = vCodeBook.get(compressedVBlocks.get(blockIndex));
                if (vVector == null) {
                    throw new IOException("Invalid V code at index " + blockIndex);
                }
                int vectorIndex = 0;
                for (int dy = 0; dy < blockSize && y + dy < height / 2; dy++) {
                    for (int dx = 0; dx < blockSize && x + dx < width / 2; dx++) {
                        vMatrix[y + dy][x + dx] = Math.min(255, Math.max(0, vVector[vectorIndex++]));
                    }
                }
                blockIndex++;
            }
        }
        
        double[][] upsampledU = upsampleChannel(uMatrix, height, width);
        double[][] upsampledV = upsampleChannel(vMatrix, height, width);
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                double Y = yMatrix[y][x];
                double U = upsampledU[y][x] - 128;
                double V = upsampledV[y][x] - 128;
                
                int r = (int) Math.round(Y + 1.140 * V);
                int g = (int) Math.round(Y - 0.395 * U - 0.581 * V);
                int b = (int) Math.round(Y + 2.032 * U);
                
                r = Math.min(255, Math.max(0, r));
                g = Math.min(255, Math.max(0, g));
                b = Math.min(255, Math.max(0, b));
                
                int rgb = (r << 16) | (g << 8) | b;
                decompressedImage.setRGB(x, y, rgb);
            }
        }
        
        return decompressedImage;
    }
    
    public static double calculateMSE(BufferedImage original, BufferedImage decompressed) {
        if (original.getWidth() != decompressed.getWidth() || original.getHeight() != decompressed.getHeight()) {
            throw new IllegalArgumentException("Images must have the same dimensions");
        }
        
        int width = original.getWidth();
        int height = original.getHeight();
        double mse = 0.0;
        long pixelCount = 0;
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int rgb1 = original.getRGB(x, y);
                int rgb2 = decompressed.getRGB(x, y);
                
                int r1 = (rgb1 >> 16) & 0xFF;
                int g1 = (rgb1 >> 8) & 0xFF;
                int b1 = rgb1 & 0xFF;
                
                int r2 = (rgb2 >> 16) & 0xFF;
                int g2 = (rgb2 >> 8) & 0xFF;
                int b2 = rgb2 & 0xFF;
                
                mse += (r1 - r2) * (r1 - r2) + (g1 - g2) * (g1 - g2) + (b1 - b2) * (b1 - b2);
                pixelCount++;
            }
        }
        
        return mse / (pixelCount * 3);
    }

    public static double calculateCompressionSize(int blockSize, int imageWidth, int imageHeight) {
        int codeBookSize = 256;
        int codeBookSize_compressed = codeBookSize * 3;
        double yBlocks = ((imageHeight * imageWidth) / (blockSize * blockSize));
        double uvBlocks = ((imageHeight / 2) * (imageWidth / 2)) / (blockSize * blockSize);
        double compressedImageSize = yBlocks + (2 * uvBlocks); 

        double total = codeBookSize_compressed + compressedImageSize;
        int originalImageSize = imageWidth * imageHeight;
        return (double) total / originalImageSize * 100;
    }

    public static void main(String[] args) {
        try {
            String folderPath = "C:\\Users\\mariam\\OneDrive\\Desktop\\vs code\\Java\\bonus-begad\\Data\\Train";
            String outputPath = "C:\\Users\\mariam\\OneDrive\\Desktop\\vs code\\Java\\bonus-begad\\Data\\Codebooks";
            String compressedOutputPath = "C:\\Users\\mariam\\OneDrive\\Desktop\\vs code\\Java\\bonus-begad\\out";
            String decompressedOutputPath = "C:\\Users\\mariam\\OneDrive\\Desktop\\vs code\\Java\\bonus-begad\\out";
            int blockSize = 2;
            int codeBookSize = 256;

            // // Create output directories if they don't exist
            // File codebookDir = new File(outputPath);
            // if (!codebookDir.exists()) {
            //     codebookDir.mkdirs();
            // }
            // File compressedDir = new File(compressedOutputPath);
            // if (!compressedDir.exists()) {
            //     compressedDir.mkdirs();
            // }
            // File decompressedDir = new File(decompressedOutputPath);
            // if (!decompressedDir.exists()) {
            //     decompressedDir.mkdirs();
            // }

            // // Generate single set of codebooks using all blocks from training images
            // List<List<int[]>> allYUVBlocks = collectAllYUVBlocks(folderPath, blockSize);
            
            // System.out.println("Total Y blocks collected: " + allYUVBlocks.get(0).size());
            // System.out.println("Total U blocks collected: " + allYUVBlocks.get(1).size());
            // System.out.println("Total V blocks collected: " + allYUVBlocks.get(2).size());

            // // Generate and save one codebook per channel
            // if (!allYUVBlocks.get(0).isEmpty()) {
            //     HashMap<String, double[]> yCodeBook = Create_Code_Book(allYUVBlocks.get(0), codeBookSize);
            //     saveCodeBook(yCodeBook, outputPath + "\\y_codebook.ser");
            //     System.out.println("Y codebook saved with " + yCodeBook.size() + " entries");
            // }
            
            // if (!allYUVBlocks.get(1).isEmpty()) {
            //     HashMap<String, double[]> uCodeBook = Create_Code_Book(allYUVBlocks.get(1), codeBookSize);
            //     saveCodeBook(uCodeBook, outputPath + "\\u_codebook.ser");
            //     System.out.println("U codebook saved with " + uCodeBook.size() + " entries");
            // }
            
            // if (!allYUVBlocks.get(2).isEmpty()) {
            //     HashMap<String, double[]> vCodeBook = Create_Code_Book(allYUVBlocks.get(2), codeBookSize);
            //     saveCodeBook(vCodeBook, outputPath + "\\v_codebook.ser");
            //     System.out.println("V codebook saved with " + vCodeBook.size() + " entries");
            // }

            // Process each image in the folder using the single set of codebooks
            String testFolderPath = "C:\\Users\\mariam\\OneDrive\\Desktop\\vs code\\Java\\bonus-begad\\Data\\Test";
            Scanner scanner = new Scanner(System.in);

            System.out.print("Enter the image file name : ");
            String fileName = scanner.nextLine();
            String filePath = testFolderPath + "\\" + fileName;

            File imageFile = new File(filePath);
            if (!imageFile.exists() || !imageFile.isFile() || 
                !(fileName.toLowerCase().endsWith(".jpg") || 
                  fileName.toLowerCase().endsWith(".png") || 
                  fileName.toLowerCase().endsWith(".jpeg"))) {
                throw new IOException("Invalid or non-existent image file: " + fileName);
            }

            System.out.println("Processing image: " + fileName);

            BufferedImage image = ImageIO.read(imageFile);

            List<List<String>> compressedChannels = Compress_test_img(image, blockSize);

            String fileNameWithoutExt = fileName.substring(0, fileName.lastIndexOf('.'));
            String yOutput = compressedOutputPath + "\\" + fileNameWithoutExt + "_y.txt";
            String uOutput = compressedOutputPath + "\\" + fileNameWithoutExt + "_u.txt";
            String vOutput = compressedOutputPath + "\\" + fileNameWithoutExt + "_v.txt";
            
            saveCompressedBlocks(compressedChannels.get(0), yOutput);
            saveCompressedBlocks(compressedChannels.get(1), uOutput);
            saveCompressedBlocks(compressedChannels.get(2), vOutput);

            BufferedImage decompressedImage = Decompress_test_img(
                yOutput,
                uOutput,
                vOutput,
                blockSize,
                image.getWidth(),
                image.getHeight()
            );

            String decompressedFile = decompressedOutputPath + "\\" + fileNameWithoutExt + "_decompressed.jpg";
            ImageIO.write(decompressedImage, "jpg", new File(decompressedFile));

            double mse = calculateMSE(image, decompressedImage);
            double compressionRatio = calculateCompressionSize(blockSize, image.getWidth(), image.getHeight());

            System.out.println("Image: " + fileName);
            System.out.println("MSE: " + mse);
            System.out.println("Compression ratio: " + compressionRatio + "%");
            System.out.println("Compressed blocks saved to: " + compressedOutputPath);
            System.out.println("Decompressed image saved to: " + decompressedFile);
            System.out.println("-----------------------------------");
            System.out.println("Image processed successfully using the existing codebooks.");

            scanner.close();

        } catch (IOException | ClassNotFoundException e) {
            e.printStackTrace();
        }
    }
}