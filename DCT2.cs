public static Matrix<double> DCT2D(Matrix<double> input)
{
    int rows = input.RowCount;
    int cols = input.ColumnCount;
    var result = Matrix<double>.Build.Dense(rows, cols);           
    // Process rows first
    var rowTemp = new Complex[cols*2]; // Double length for symmetry
    for (int i = 0; i < rows; i++)
    {
        // Create symmetric sequence
        for (int j = 0; j < cols; j++)
        {
            rowTemp[j] = new Complex(input[i, j], 0);
            rowTemp[cols*2 - 1 - j] = new Complex(input[i, j], 0);
        }
       
        // Compute FFT
        Fourier.Forward(rowTemp, FourierOptions.Matlab);               

        // Extract DCT coefficients
        double scale = Math.Sqrt(2.0 * cols);
        for (int k = 0; k < cols; k++)
        {
            double phase = -Math.PI * k / (2.0 * cols);
            var c = Complex.Exp(new Complex(0, phase));
            result[i, k] = (rowTemp[k] * c).Real / scale * (k == 0 ? 1.0 / Math.Sqrt(2) : 1.0);
        }
    }
   
    var colTemp = new Complex[rows*2];
    for (int j = 0; j < cols; j++)
    {
        // Create symmetric sequence
        for (int i = 0; i < rows; i++)
        {
            colTemp[i] = new Complex(result[i, j], 0);
            colTemp[rows*2 - 1 - i] = new Complex(result[i, j], 0);                   
        }

        // Compute FFT
        Fourier.Forward(colTemp, FourierOptions.Matlab);

        // Extract DCT coefficients
        double scale = Math.Sqrt(2.0 * rows);
        for (int k = 0; k < rows; k++)
        {
            double phase = -Math.PI * k / (2.0 * rows);
            var c = Complex.Exp(new Complex(0, phase));
            result[k, j] = (colTemp[k] * c).Real / scale * (k == 0 ? 1.0 / Math.Sqrt(2) : 1.0);
        }
    }

    return result;
}