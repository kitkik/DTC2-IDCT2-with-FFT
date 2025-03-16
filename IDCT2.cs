 public static Matrix<double> IDCT2D(Matrix<double> input)
 {
     int rows = input.RowCount;
     int cols = input.ColumnCount;
     var result = Matrix<double>.Build.Dense(rows, cols);            

     // Process columns first
     var colTemp = new Complex[2 * rows];
     for (int j = 0; j < cols; j++)
     {
         // Prepare input with proper scaling
         for (int k = 0; k < rows; k++)
         {
             double scale = Math.Sqrt(2.0 / rows) * (k == 0 ? 1.0/Math.Sqrt(2) : 1.0);
             double phase = Math.PI * k / (2.0 * rows);
             var c = Complex.Exp(new Complex(0, phase));
             colTemp[k] = (input[k, j] * scale * c);
             colTemp[2 * rows - 1 - k] = new Complex(0, 0); // Will be filled by IFFT
         }
         // Compute IFFT
         Fourier.Inverse(colTemp, FourierOptions.Matlab);

         // Extract real part
         for (int i = 0; i < rows; i++)
         {
             result[i, j] = colTemp[i].Real * 2 * rows;
         }
     }

     // Process rows
     var rowTemp = new Complex[2 * cols];
     for (int i = 0; i < rows; i++)
     {
         // Prepare input with proper scaling
         for (int k = 0; k < cols; k++)
         {
             double scale = Math.Sqrt(2.0 / cols) * (k == 0 ? 1.0 / Math.Sqrt(2) : 1.0);
             double phase = Math.PI * k / (2.0 * cols);
             var c = Complex.Exp(new Complex(0, phase));
             rowTemp[k] = (result[i, k] * scale* c);
             rowTemp[2 * cols - 1 - k] = new Complex(0, 0);
         }
         // Compute IFFT
         Fourier.Inverse(rowTemp, FourierOptions.Matlab);

         // Extract real part
         for (int j = 0; j < cols; j++)
         {
             result[i, j] = rowTemp[j].Real * 2 * cols;
         }
     }

     return result;
 }