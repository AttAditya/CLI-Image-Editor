import java.util.*;

import java.io.File;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.awt.Color;

import java.io.IOException;

import java.math.*;

class ImageEditor {
    static int limitColor0To255(double valueRaw) {
        return (int) Math.max(0, Math.min(255, valueRaw));
    }

    static int[][] imageToMatrix(BufferedImage buffImage) {
        int width = buffImage.getWidth();
        int height = buffImage.getHeight();

        int[][] matrix = new int[height][width];

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                matrix[y][x] = buffImage.getRGB(x, y);
            }
        }

        return matrix;
    }

    static int matrixSum(int[][] mtrx) {
        int sum = 0;

        for (int[] arr: mtrx) {
            for (int el: arr) {
                sum += el;
            }
        }

        return sum;
    }

    static double matrixMean(int[][] mtrx) {
        double sum = matrixSum(mtrx);
        int count = mtrx.length * mtrx[0].length;

        return sum / count;
    }

    static void swapMatrixElements(int[][] mtrx, int i1, int j1, int i2, int j2) {
        int temp = mtrx[i1][j1];
        mtrx[i1][j1] = mtrx[i2][j2];
        mtrx[i2][j2] = temp;
    }

    static void flipMatrixHorizontally(int[][] mtrx) {
        for (int i = 0; i < mtrx.length; i++) {
            for (int j = 0; j < mtrx[0].length / 2; j++) {
                swapMatrixElements(mtrx, i, j, i, mtrx[0].length - j - 1);
            }
        }
    }

    static void flipMatrixVertically(int[][] mtrx) {
        for (int i = 0; i < mtrx.length / 2; i++) {
            for (int j = 0; j < mtrx[0].length; j++) {
                swapMatrixElements(mtrx, i, j, mtrx.length - i - 1, j);
            }
        }
    }

    static int[][] transposeMatrix(int[][] matrix) {
        int[][] mtrx = new int[matrix[0].length][matrix.length];

        for (int i = 0; i < mtrx.length; i++) {
            for (int j = 0; j < mtrx[0].length; j++) {
                mtrx[i][j] = matrix[j][i];
            }
        }

        return mtrx;
    }

    static int[][] matrixCircle(int[][] matrix) {
        int radius = Math.min(matrix.length, matrix[0].length) / 2;
        int[][] resultMatrix = new int[matrix.length][matrix[0].length];

        int xOffset = (matrix[0].length / 2) - radius;
        int yOffset = (matrix.length / 2) - radius;

        for (int x = 0; x < matrix[0].length; x++) {
            for (int y = 0; y < matrix.length; y++) {
                if (((x - radius - xOffset) * (x - radius - xOffset)) + ((y - radius - yOffset) * (y - radius - yOffset)) <= (radius * radius)) {
                    resultMatrix[y][x] = matrix[y][x];
                }
            }
        }

        return resultMatrix;
    }

    static int[][] matrixVignette(int[][] matrix) {
        int radius = Math.min(matrix.length, matrix[0].length) / 2;
        int[][] resultMatrix = new int[matrix.length][matrix[0].length];

        int xOffset = (matrix[0].length / 2) - radius;
        int yOffset = (matrix.length / 2) - radius;

        for (int x = 0; x < matrix[0].length; x++) {
            for (int y = 0; y < matrix.length; y++) {
                double xValSq = ((x - radius - xOffset) * (x - radius - xOffset));
                double yValSq = ((y - radius - yOffset) * (y - radius - yOffset));
                double rValSq = (radius * radius);

                Color currentColor = new Color(matrix[y][x]);

                double distCenter = xValSq + yValSq;
                double distCirleBound = distCenter - rValSq;

                if (distCirleBound < 0) {
                    distCirleBound = 1;
                }

                int rVal = (int) (currentColor.getRed() / Math.sqrt(distCirleBound));
                int gVal = (int) (currentColor.getGreen() / Math.sqrt(distCirleBound));
                int bVal = (int) (currentColor.getBlue() / Math.sqrt(distCirleBound));

                rVal = limitColor0To255(rVal);
                gVal = limitColor0To255(gVal);
                bVal = limitColor0To255(bVal);

                Color updatedColor = new Color(rVal, gVal, bVal);

                resultMatrix[y][x] = updatedColor.getRGB();
            }
        }

        return resultMatrix;
    }

    static int[][] subMatrix(int[][] matrix, int x1, int y1, int x2, int y2) {
        int[][] resultMatrix = new int[y2 - y1 + 1][x2 - x1 + 1];

        for (int x = x1; x <= x2; x++) {
            for (int y = y1; y <= y2; y++) {
                if (x >= matrix[0].length || y >= matrix.length || x < 0 || y < 0) {
                    continue;
                }

                resultMatrix[y - y1][x - x1] = matrix[y][x];
            }
        }

        return resultMatrix;
    }

    static int[][][] splitRGB(int[][] matrix) {
        int[][][] splitMatrix = new int[3][matrix.length][matrix[0].length];

        for (int y = 0; y < matrix.length; y++) {
            for (int x = 0; x < matrix[0].length; x++) {
                Color colors = new Color(matrix[y][x]);
                splitMatrix[0][y][x] = colors.getRed();
            }
        }

        for (int y = 0; y < matrix.length; y++) {
            for (int x = 0; x < matrix[0].length; x++) {
                Color colors = new Color(matrix[y][x]);
                splitMatrix[1][y][x] = colors.getGreen();
            }
        }

        for (int y = 0; y < matrix.length; y++) {
            for (int x = 0; x < matrix[0].length; x++) {
                Color colors = new Color(matrix[y][x]);
                splitMatrix[2][y][x] = colors.getBlue();
            }
        }

        return splitMatrix;
    }

    static int[][] contrastMatrix(int[][] matrix, int contrastAmount) {
        int[][] resultMatrix = new int[matrix.length][matrix[0].length];

        for (int x = 0; x < matrix[0].length; x++) {
            for (int y = 0; y < matrix.length; y++) {
                Color currentColor = new Color(matrix[y][x]);

                int rValue = currentColor.getRed();
                int gValue = currentColor.getGreen();
                int bValue = currentColor.getBlue();

                if (rValue > 127) {
                    rValue += contrastAmount;
                } else {
                    rValue -= contrastAmount;
                }
                
                if (gValue > 127) {
                    gValue += contrastAmount;
                } else {
                    gValue -= contrastAmount;
                }

                if (bValue > 127) {
                    bValue += contrastAmount;
                } else {
                    bValue -= contrastAmount;
                }

                Color updatedColor = new Color(limitColor0To255(rValue), limitColor0To255(gValue), limitColor0To255(bValue));

                resultMatrix[y][x] = updatedColor.getRGB();
            }
        }

        return resultMatrix;
    }

    static int[][] colorInverseMatrix(int[][] matrix) {
        int[][] resultMatrix = new int[matrix.length][matrix[0].length];

        for (int x = 0; x < matrix[0].length; x++) {
            for (int y = 0; y < matrix.length; y++) {
                Color currentColor = new Color(matrix[y][x]);

                int rValue = 255 - currentColor.getRed();
                int gValue = 255 - currentColor.getGreen();
                int bValue = 255 - currentColor.getBlue();

                Color updatedColor = new Color(limitColor0To255(rValue), limitColor0To255(gValue), limitColor0To255(bValue));

                resultMatrix[y][x] = updatedColor.getRGB();
            }
        }

        return resultMatrix;
    }

    static int[] rotateCoordinate(int x, int y, double sinAngle, double cosAngle) {
        int updatedX = (int) ((x * cosAngle) - (y * sinAngle));
        int updatedY = (int) ((y * cosAngle) + (x * sinAngle));

        int[] coords = {updatedX, updatedY};

        return coords;
    }

    static int[][] rotateMatrix(int[][] matrix, int angle) {
        angle = angle % 360;

        if (angle < 0) {
            angle += 360;
        }

        double sinAngle = Math.sin(angle * Math.PI / 180);
        double cosAngle = Math.cos(angle * Math.PI / 180);

        int cornerTLX = rotateCoordinate(0, 0, sinAngle, cosAngle)[0];
        int cornerTLY = rotateCoordinate(0, 0, sinAngle, cosAngle)[1];

        int cornerTRX = rotateCoordinate(matrix[0].length - 1, 0, sinAngle, cosAngle)[0];
        int cornerTRY = rotateCoordinate(matrix[0].length - 1, 0, sinAngle, cosAngle)[1];

        int cornerBLX = rotateCoordinate(0, matrix.length - 1, sinAngle, cosAngle)[0];
        int cornerBLY = rotateCoordinate(0, matrix.length - 1, sinAngle, cosAngle)[1];

        int cornerBRX = rotateCoordinate(matrix[0].length - 1, matrix.length - 1, sinAngle, cosAngle)[0];
        int cornerBRY = rotateCoordinate(matrix[0].length - 1, matrix.length - 1, sinAngle, cosAngle)[1];

        int[] xCoords = {cornerBLX, cornerTLX, cornerBRX, cornerTRX};
        int[] yCoords = {cornerBLY, cornerTLY, cornerBRY, cornerTRY};

        Arrays.sort(xCoords);
        Arrays.sort(yCoords);

        int xOffset = -xCoords[0];
        int yOffset = -yCoords[0];

        int width = xCoords[3] - xCoords[0];
        int height = yCoords[3] - yCoords[0];

        int[][] resultMatrix = new int[height][width];

        for (int x = 0; x < matrix[0].length; x++) {
            for (int y = 0; y < matrix.length; y++) {
                int updatedX = rotateCoordinate(x, y, sinAngle, cosAngle)[0];
                int updatedY = rotateCoordinate(x, y, sinAngle, cosAngle)[1];

                updatedX += xOffset;
                updatedY += yOffset;

                if (updatedX < 0 || updatedX >= resultMatrix[0].length || updatedY < 0 || updatedY >= resultMatrix.length) {
                    continue;
                }

                resultMatrix[updatedY][updatedX] = matrix[y][x];
            }
        }

        return resultMatrix;
    }

    static int[][] binaryColorMatrix(int[][] matrix, double tolerance) {
        int[][] resultMatrix = new int[matrix.length][matrix[0].length];

        for (int x = 0; x < matrix[0].length; x++) {
            for (int y = 0; y < matrix.length; y++) {
                Color currentColor = new Color(matrix[y][x]);

                int rValue = currentColor.getRed();
                int gValue = currentColor.getGreen();
                int bValue = currentColor.getBlue();

                if (rValue > 255 * (1 - tolerance)) {
                    rValue = 255;
                } else {
                    rValue = 0;
                }
                
                if (gValue > 255 * (1 - tolerance)) {
                    gValue = 255;
                } else {
                    gValue = 0;
                }

                if (bValue > 255 * (1 - tolerance)) {
                    bValue = 255;
                } else {
                    bValue = 0;
                }

                Color updatedColor = new Color(limitColor0To255(rValue), limitColor0To255(gValue), limitColor0To255(bValue));

                resultMatrix[y][x] = updatedColor.getRGB();
            }
        }

        return resultMatrix;
    }

    static int[][] blurMatrix(int[][] matrix, int blurRadius) {
        int[][] resultMatrix = new int[matrix.length][matrix[0].length];

        for (int x = 0; x < matrix[0].length; x++) {
            for (int y = 0; y < matrix.length; y++) {
                int[][] subMtrx = subMatrix(matrix, x - (blurRadius / 2), y - (blurRadius / 2), x + (blurRadius / 2), y + (blurRadius / 2));

                int[][][] RGBMatrices = splitRGB(subMtrx);

                int[][] RMatrix = RGBMatrices[0];
                int[][] GMatrix = RGBMatrices[1];
                int[][] BMatrix = RGBMatrices[2];

                int RMean = (int) matrixMean(RMatrix);
                int GMean = (int) matrixMean(GMatrix);
                int BMean = (int) matrixMean(BMatrix);

                Color updatedColor = new Color(RMean, GMean, BMean);

                resultMatrix[y][x] = updatedColor.getRGB();
            }
        }

        return resultMatrix;
    }

    static int[][] pixelateMatrix(int[][] matrix, int blurRadius) {
        int[][] resultMatrix = new int[matrix.length][matrix[0].length];

        for (int x = 0; x < matrix[0].length; x++) {
            for (int y = 0; y < matrix.length; y++) {
                int[][] subMtrx = subMatrix(matrix, (x - (x % blurRadius)) - (blurRadius / 2), (y - (y % blurRadius)) - (blurRadius / 2), (x - (x % blurRadius)) + (blurRadius / 2), (y - (y % blurRadius)) + (blurRadius / 2));

                int[][][] RGBMatrices = splitRGB(subMtrx);

                int[][] RMatrix = RGBMatrices[0];
                int[][] GMatrix = RGBMatrices[1];
                int[][] BMatrix = RGBMatrices[2];

                int RMean = (int) matrixMean(RMatrix);
                int GMean = (int) matrixMean(GMatrix);
                int BMean = (int) matrixMean(BMatrix);

                Color updatedColor = new Color(RMean, GMean, BMean);

                resultMatrix[y][x] = updatedColor.getRGB();
            }
        }

        return resultMatrix;
    }

    static int[][] edgeDetectMatrix(int[][] matrix, int edgeWidth) {
        int[][] resultMatrix = new int[matrix.length][matrix[0].length];

        for (int x = 0; x < matrix[0].length - edgeWidth; x++) {
            for (int y = 0; y < matrix.length - edgeWidth; y++) {
                Color currentColor = new Color(matrix[y][x]);
                Color ShiftedColor = new Color(matrix[y + edgeWidth][x + edgeWidth]);

                int rValue = currentColor.getRed() - ShiftedColor.getRed();
                int gValue = currentColor.getGreen() - ShiftedColor.getGreen();
                int bValue = currentColor.getBlue() - ShiftedColor.getBlue();

                Color updatedColor = new Color(limitColor0To255(rValue), limitColor0To255(gValue), limitColor0To255(bValue));

                resultMatrix[y][x] = updatedColor.getRGB();
            }
        }

        return resultMatrix;
    }

    static int[][] edgeDetectMatrixV2(int[][] matrix, int edgeWidth) {
        int[][] resultMatrix = new int[matrix.length][matrix[0].length];
        int[][] resultMatrix1 = new int[matrix.length][matrix[0].length];
        int[][] resultMatrix2 = new int[matrix.length][matrix[0].length];

        for (int x = 0; x < matrix[0].length - edgeWidth; x++) {
            for (int y = 0; y < matrix.length - edgeWidth; y++) {
                Color currentColor = new Color(matrix[y][x]);
                Color ShiftedColor = new Color(matrix[y + edgeWidth][x + edgeWidth]);

                int rValue = currentColor.getRed() - ShiftedColor.getRed();
                int gValue = currentColor.getGreen() - ShiftedColor.getGreen();
                int bValue = currentColor.getBlue() - ShiftedColor.getBlue();

                Color updatedColor = new Color(limitColor0To255(rValue), limitColor0To255(gValue), limitColor0To255(bValue));

                resultMatrix1[y][x] = updatedColor.getRGB();
            }
        }

        for (int x = edgeWidth; x < matrix[0].length; x++) {
            for (int y = edgeWidth; y < matrix.length - edgeWidth; y++) {
                Color currentColor = new Color(matrix[y][x]);
                Color ShiftedColor = new Color(matrix[y - edgeWidth][x - edgeWidth]);

                int rValue = currentColor.getRed() - ShiftedColor.getRed();
                int gValue = currentColor.getGreen() - ShiftedColor.getGreen();
                int bValue = currentColor.getBlue() - ShiftedColor.getBlue();

                Color updatedColor = new Color(limitColor0To255(rValue), limitColor0To255(gValue), limitColor0To255(bValue));

                resultMatrix2[y][x] = updatedColor.getRGB();
            }
        }

        for (int x = 0; x < matrix[0].length; x++) {
            for (int y = 0; y < matrix.length; y++) {
                Color matrix1Color = new Color(resultMatrix1[y][x]);
                Color matrix2Color = new Color(resultMatrix2[y][x]);

                int rValue = matrix1Color.getRed() + matrix2Color.getRed();
                int gValue = matrix1Color.getGreen() + matrix2Color.getGreen();
                int bValue = matrix1Color.getBlue() + matrix2Color.getBlue();

                Color updatedColor = new Color(limitColor0To255(rValue), limitColor0To255(gValue), limitColor0To255(bValue));

                resultMatrix[y][x] = updatedColor.getRGB();
            }
        }

        return resultMatrix;
    }

    static BufferedImage matrixToBufferedImage(int[][] matrix) {
        int width = matrix[0].length;
        int height = matrix.length;

        BufferedImage outImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                Color rgbValues = new Color(matrix[y][x]);
                Color newColors = new Color(rgbValues.getRed(), rgbValues.getGreen(), rgbValues.getBlue());
                
                outImage.setRGB(x, y, newColors.getRGB());
            }
        }

        return outImage;
    }

    static BufferedImage grayscale(BufferedImage buffImage) {
        int width = buffImage.getWidth();
        int height = buffImage.getHeight();

        BufferedImage outImage = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                Color rgbValues = new Color(buffImage.getRGB(x, y));
                int grayValue = (rgbValues.getRed() + rgbValues.getGreen() + rgbValues.getBlue()) / 3;
                outImage.setRGB(x, y, rgbValues.getRGB());
            }
        }

        return outImage;
    }

    static BufferedImage colorFixImage(BufferedImage buffImage, double redRatio, double greenRatio, double blueRatio) {
        int width = buffImage.getWidth();
        int height = buffImage.getHeight();

        BufferedImage outImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                Color rgbValues = new Color(buffImage.getRGB(x, y));
                Color newColors = new Color(limitColor0To255(rgbValues.getRed() * redRatio), limitColor0To255(rgbValues.getGreen() * greenRatio), limitColor0To255(rgbValues.getBlue() * blueRatio));
                outImage.setRGB(x, y, newColors.getRGB());
            }
        }

        return outImage;
    }

    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);

        while (true) {
            try {
                System.out.println();
                System.out.print("Filename(Leave empty to exit): ");
                String filename = sc.nextLine();
   
                System.out.println();

                if (filename.length() == 0) {
                    return;
                }

                File inputFile = new File(filename);
                BufferedImage buffImage = ImageIO.read(inputFile);
   
                System.out.println("(0)\tCancel");
                System.out.println("(1)\tGrayscale");
                System.out.println("(2)\tRotate Clockwise");
                System.out.println("(3)\tRotate Anti-Clockwise");
                System.out.println("(4)\tFlip Horizontally");
                System.out.println("(5)\tFlip Vertically");
                System.out.println("(6)\tBrightness");
                System.out.println("(7)\tBlur");
                System.out.println("(8)\tPixelate");
                System.out.println("(9)\tCrop to circle");
                System.out.println("(10)\tContrast");
                System.out.println("(11)\tColor Inversion");
                System.out.println("(12)\tEdge Detection");
                System.out.println("(13)\tSketch");
                System.out.println("(14)\tRotate");

                System.out.println();
                System.out.print("Option: ");
   
                int action = sc.nextInt();
                sc.nextLine();

                if (action == 0) {
                    System.out.println("Cancelled operation on image...");
                    continue;
                }

                System.out.print("Output filename: ");
   
                String outputFilename = sc.nextLine();

                int[][] imageMatrix = imageToMatrix(buffImage);

                File outputFile = new File(outputFilename);
                double amount;

                switch (action) {
                    case 1:
                        ImageIO.write(grayscale(buffImage), "jpg", outputFile);
   
                        break;
                    case 2:
                        imageMatrix = transposeMatrix(imageMatrix);
                        flipMatrixHorizontally(imageMatrix);

                        ImageIO.write(matrixToBufferedImage(imageMatrix), "jpg", outputFile);

                        break;
                    case 3:
                        for (int i = 1; i <= 3; i++) {
                            imageMatrix = transposeMatrix(imageMatrix);
                            flipMatrixHorizontally(imageMatrix);
                        }

                        ImageIO.write(matrixToBufferedImage(imageMatrix), "jpg", outputFile);

                        break;
                    case 4:
                        flipMatrixHorizontally(imageMatrix);
                        ImageIO.write(matrixToBufferedImage(imageMatrix), "jpg", outputFile);

                        break;
                    case 5:
                        flipMatrixVertically(imageMatrix);
                        ImageIO.write(matrixToBufferedImage(imageMatrix), "jpg", outputFile);

                        break;
                    case 6:
                        System.out.print("Amount(in percentage / 100 range, ex - 120% -> 1.2, 80% -> 0.8): ");
                        amount = sc.nextDouble();
                        sc.nextLine();
                        
                        ImageIO.write(colorFixImage(buffImage, amount, amount, amount), "jpg", outputFile);

                        break;
                    case 7:
                        imageMatrix = blurMatrix(imageMatrix, 10);
                        ImageIO.write(matrixToBufferedImage(imageMatrix), "jpg", outputFile);

                        break;
                    case 8:
                        imageMatrix = pixelateMatrix(imageMatrix, 10);
                        ImageIO.write(matrixToBufferedImage(imageMatrix), "jpg", outputFile);

                        break;
                    case 9:
                        imageMatrix = matrixCircle(imageMatrix);
                        ImageIO.write(matrixToBufferedImage(imageMatrix), "jpg", outputFile);

                        break;
                    case 10:
                        System.out.print("Amount(in 0.0 -> 1.0 range): ");
                        amount = sc.nextDouble();
                        sc.nextLine();

                        imageMatrix = contrastMatrix(imageMatrix, (int) (255 * amount));
                        ImageIO.write(matrixToBufferedImage(imageMatrix), "jpg", outputFile);

                        break;
                    case 11:
                        imageMatrix = colorInverseMatrix(imageMatrix);
                        ImageIO.write(matrixToBufferedImage(imageMatrix), "jpg", outputFile);

                        break;
                    case 12:
                        imageMatrix = edgeDetectMatrix(imageMatrix, 5);
                        ImageIO.write(matrixToBufferedImage(imageMatrix), "jpg", outputFile);

                        break;
                    case 13:
                        System.out.print("Stroke Amount(in pixels, recommended 2px): ");
                        int stroke = sc.nextInt();
                        sc.nextLine();

                        System.out.print("Tolerance Amount(in 0.0 -> 1.0 range, recommended 0.97): ");
                        double tolerance = sc.nextDouble();
                        sc.nextLine();

                        buffImage = grayscale(buffImage);
                        
                        imageMatrix = imageToMatrix(buffImage);
                        imageMatrix = edgeDetectMatrix(imageMatrix, stroke);
                        imageMatrix = binaryColorMatrix(imageMatrix, tolerance);
                        imageMatrix = colorInverseMatrix(imageMatrix);

                        ImageIO.write(matrixToBufferedImage(imageMatrix), "jpg", outputFile);

                        break;
                    case 14:
                        System.out.print("Angle(in degrees): ");
                        int angle = sc.nextInt();
                        sc.nextLine();

                        imageMatrix = rotateMatrix(imageMatrix, angle);
                        ImageIO.write(matrixToBufferedImage(imageMatrix), "jpg", outputFile);

                        break;
                    default:
                        System.out.println("Option invalid!");
                        
                        break;
                }
            } catch (IOException err) {
                System.out.println("Please enter a valid file...");
            }
        }
    }
}
