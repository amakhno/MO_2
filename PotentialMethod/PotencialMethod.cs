using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PotentialMethod
{
    class PotencialMethod
    {
        Matrix<double> PriceMatrixs;
        Matrix<double> CurrentStepMatrix;
        static MatrixBuilder<double> M = Matrix<double>.Build;
        enum Status { Empty = -1, Negative = -2, Positive = -3, NotInit = -4 }

        public PotencialMethod(double[,] inputMatrix, Matrix<double> priviousStepMatrix)
        {
            this.CurrentStepMatrix = priviousStepMatrix;
            this.PriceMatrixs = M.DenseOfArray(inputMatrix);

            string alphaindexes = "";
            string betaindexes = "";

            while (true)
            {
                for (int i = 0; i < priviousStepMatrix.RowCount; i++)
                {
                    for (int j = 0; j < priviousStepMatrix.ColumnCount; j++)
                    {
                        if (priviousStepMatrix[i, j] != (double)Status.Empty)
                        {
                            alphaindexes += i.ToString();
                            betaindexes += j.ToString();
                        }
                    }
                }
                var alpha = alphaindexes.OrderBy(x => x.ToString()).Distinct();
                var beta = betaindexes.OrderBy(x => x.ToString()).Distinct();
                alphaindexes = "";
                betaindexes = "";
                foreach (char a in alpha)
                {
                    alphaindexes += a;
                }
                foreach (char a in beta)
                {
                    betaindexes += a;
                }
                var Matrix = M.Sparse(alphaindexes.Length + betaindexes.Length, alphaindexes.Length + betaindexes.Length);
                double[] rightSide = new double[Matrix.ColumnCount];
                rightSide[0] = 0;
                Matrix[0, 0] = 1;
                int currentRow = 1;
                for (int i = 0; i < priviousStepMatrix.RowCount; i++)
                {
                    for (int j = 0; j < priviousStepMatrix.ColumnCount; j++)
                    {
                        if (priviousStepMatrix[i, j] != (double)Status.Empty)
                        {
                            Matrix[currentRow, i] = 1;
                            Matrix[currentRow, j + alphaindexes.Length] = 1;
                            rightSide[currentRow] = PriceMatrixs[i, j];
                            currentRow++;
                        }
                    }
                }
                var solve = Matrix.Solve(DenseVector.Build.DenseOfArray(rightSide)).ToArray();
                double currentMin = double.MaxValue;
                int[] AdressMin = new int[] { 0 };
                for (int i = 0; i < priviousStepMatrix.RowCount; i++)
                {
                    for (int j = 0; j < priviousStepMatrix.ColumnCount; j++)
                    {
                        if (priviousStepMatrix[i, j] == (double)Status.Empty)
                        {
                            if ((PriceMatrixs[i, j] - (i + j)) < currentMin)
                            {
                                currentMin = (PriceMatrixs[i, j] - (i + j));
                                AdressMin = new int[] { i, j };
                            }
                        }
                    }
                }
                //Console.Write(priviousStepMatrix.ToString());
                var Cycle = BuildCycle(priviousStepMatrix, AdressMin[0], AdressMin[1]);
            }
        }

        private List<int> BuildCycle(Matrix<double> priviousStepMatrix, int i, int j)
        {
            List<int> way = new List<int>();
            way.Add(i);
            way.Add(j);
            priviousStepMatrix[i, j] = (double)Status.Positive;
            var a = GoRow(priviousStepMatrix, i, j, way, false);
            if (a.Count() != 0) return a;
            a = GoColumn(priviousStepMatrix, i, j, way, false);
            if (a.Count() != 0)
                return a;
            return null;
        }

        private List<int> GoRow(Matrix<double> priviousStepMatrixOld, int i, int j, List<int> currentWayOld, bool IsPositive)
        {
            var priviousStepMatrix = M.Dense(priviousStepMatrixOld.RowCount, priviousStepMatrixOld.ColumnCount);
            priviousStepMatrixOld.CopyTo(priviousStepMatrix);
            var currentWay = currentWayOld;
            //Console.Write(priviousStepMatrix.ToString());
            for (int k = 0; k < priviousStepMatrix.RowCount; k++)
            {
                if (i == k)
                {
                    continue;
                }
                else
                {
                    if ((priviousStepMatrix[k, j] != (double)Status.Empty) && (priviousStepMatrix[k, j] != (double)Status.Positive) && (priviousStepMatrix[k, j] != (double)Status.Negative))
                    {
                        currentWay.Add(k);
                        currentWay.Add(j);
                        if (IsPositive)
                        {
                            priviousStepMatrix[k, j] = (double)Status.Positive;
                        }
                        else
                        {
                            priviousStepMatrix[k, j] = (double)Status.Negative;
                        }
                        var a = GoColumn(priviousStepMatrix, k, j, currentWay, !IsPositive);
                        if (a != null)
                            return a;
                    }
                }
            }
            if ((i == currentWay[0]) && (j == currentWay[1]))
                return currentWay;
            return null;
        }

        private List<int> GoColumn(Matrix<double> priviousStepMatrixOld, int i, int j, List<int> currentWayOld, bool IsPositive)
        {
            var priviousStepMatrix = M.Dense(priviousStepMatrixOld.RowCount, priviousStepMatrixOld.ColumnCount);
            priviousStepMatrixOld.CopyTo(priviousStepMatrix);
            var currentWay = currentWayOld;
            //Console.Write(priviousStepMatrix.ToString());
            for (int k = 0; k < priviousStepMatrix.ColumnCount; k++)
            {
                if (i == k)
                {
                    continue;
                }
                else
                {
                    if ((priviousStepMatrix[i, k] != (double)Status.Empty) && (priviousStepMatrix[i, k] != (double)Status.Positive) && (priviousStepMatrix[i, k] != (double)Status.Negative))
                    {
                        currentWay.Add(i);
                        currentWay.Add(k);
                        if (IsPositive)
                        {
                            priviousStepMatrix[i, k] = (double)Status.Positive;
                        }
                        else
                        {
                            priviousStepMatrix[i, k] = (double)Status.Negative;
                        }
                        if (GoRow(priviousStepMatrix, i, k, currentWay, !IsPositive) != null)
                            return GoRow(priviousStepMatrix, i, k, currentWay, !IsPositive);
                    }
                }
            }
            if ((i == currentWay[0]) && (j == currentWay[1]))
                return currentWay;
            return null;
        }

        static public Matrix<double> BiuldFirstStep(double[,] inputArray, double[] CanGiveArray/*N*/, double[] NeedArray/*M*/)
        {
            Matrix<double> PriviousStepMatrix = M.Dense(inputArray.GetLength(0), inputArray.GetLength(1), (double)Status.NotInit);
            int target = inputArray.GetLength(0) + inputArray.GetLength(1) - 1 - 2;
            Matrix<double> PriceTempMatrix = M.DenseOfArray(inputArray);
            int step = 0;
            while ((CanGiveArray.Sum() + NeedArray.Sum()) != 0)
            {
                step++;

                int[] AdressMin = Min(ref PriceTempMatrix, PriviousStepMatrix);
                /*if (step == target)
                {
                    AdressMin = new int[] { Max(CanGiveArray), Max(CanGiveArray) };
                }*/
                PriceTempMatrix[AdressMin[0], AdressMin[1]] = (double)Status.Empty;
                if (CanGiveArray[AdressMin[0]] >= NeedArray[AdressMin[1]])
                {
                    PriviousStepMatrix[AdressMin[0], AdressMin[1]] = NeedArray[AdressMin[1]];
                    for (int i = 0; i < CanGiveArray.Length; i++)
                    {
                        if (i == AdressMin[0])
                        {
                            continue;
                        }
                        if (PriviousStepMatrix[i, AdressMin[1]] != (double)Status.NotInit)
                            continue;
                        PriviousStepMatrix[i, AdressMin[1]] = (double)Status.Empty;
                    }
                }
                else
                {
                    PriviousStepMatrix[AdressMin[0], AdressMin[1]] = CanGiveArray[AdressMin[0]];
                    for (int i = 0; i < NeedArray.Length; i++)
                    {
                        if (i == AdressMin[1])
                        {
                            continue;
                        }
                        if (PriviousStepMatrix[AdressMin[0], i] != (double)Status.NotInit)
                            continue;
                        PriviousStepMatrix[AdressMin[0], i] = (double)Status.Empty;
                    }
                }
                CanGiveArray[AdressMin[0]] -= PriviousStepMatrix[AdressMin[0], AdressMin[1]];
                NeedArray[AdressMin[1]] -= PriviousStepMatrix[AdressMin[0], AdressMin[1]];
            }
            //Console.WriteLine("Матрица поставок");
            //Console.Write(PriviousStepMatrix.ToString());
            return PriviousStepMatrix;
        }

        static private int[] Max(ref Matrix<double> PriceTempMatrix)
        {
            int maxI = -1, maxJ = -1;
            double MaxElement = double.MinValue;
            for (int i = 0; i < PriceTempMatrix.RowCount; i++)
            {
                for (int j = 0; j < PriceTempMatrix.ColumnCount; j++)
                {
                    if (PriceTempMatrix[i, j] > MaxElement)
                    {
                        maxI = i;
                        maxJ = j;
                        MaxElement = PriceTempMatrix[i, j];
                    }
                }
            }
            return new int[] { maxI, maxJ };
        }

        static private int[] Min(ref Matrix<double> PriceTempMatrix, Matrix<double> PriviousStepMatrix)
        {
            int minI = -1, minJ = -1;
            double MinElement = double.MaxValue;
            for (int i = 0; i < PriceTempMatrix.RowCount; i++)
            {
                for (int j = 0; j < PriceTempMatrix.ColumnCount; j++)
                {
                    if ((PriceTempMatrix[i, j] < MinElement) && (PriviousStepMatrix[i, j] != (double)Status.Empty) && (PriceTempMatrix[i, j] != (double)Status.Empty))
                    {
                        minI = i;
                        minJ = j;
                        MinElement = PriceTempMatrix[i, j];
                    }
                }
            }
            return new int[] { minI, minJ };
        }

        static private double Max(double a, double b)
        {
            if (a > b) return a;
            else return b;
        }

        static private int Max(double[] array)
        {
            for (int i = 0; i < array.Length; i++)
            {
                if (array[i] == array.Min()) return i;
            }
            return -1;
        }

        static public Matrix<double> BiuldFirstStep2(double[,] inputArray, double[] CanGiveArray/*N*/, double[] NeedArray/*M*/)
        {
            Matrix<double> PriviousStepMatrix = M.Dense(inputArray.GetLength(0), inputArray.GetLength(1), (double)Status.NotInit);
            int target = inputArray.GetLength(0) + inputArray.GetLength(1) - 1 - 2;
            Matrix<double> PriceTempMatrix = M.DenseOfArray(inputArray);
            int step = 0;
            for (int i = 0; i < PriviousStepMatrix.RowCount;)
            {
                for (int j = 0; j < PriviousStepMatrix.ColumnCount;)
                {
                    int[] AdressMin = new int[] { i, j };
                    /*if (step == target)
                    {
                        AdressMin = new int[] { Max(CanGiveArray), Max(CanGiveArray) };
                    }*/
                    PriceTempMatrix[AdressMin[0], AdressMin[1]] = (double)Status.Empty;
                    if (CanGiveArray[AdressMin[0]] >= NeedArray[AdressMin[1]])
                    {
                        PriviousStepMatrix[AdressMin[0], AdressMin[1]] = NeedArray[AdressMin[1]];
                        for (int k = 0; k < CanGiveArray.Length; k++)
                        {
                            if (k == AdressMin[0])
                            {
                                continue;
                            }
                            if (PriviousStepMatrix[k, AdressMin[1]] != (double)Status.NotInit)
                                continue;
                            PriviousStepMatrix[k, AdressMin[1]] = (double)Status.Empty;
                        }
                        j++;
                        step++;
                    }
                    else
                    {
                        PriviousStepMatrix[AdressMin[0], AdressMin[1]] = CanGiveArray[AdressMin[0]];
                        for (int k = 0; k < NeedArray.Length; k++)
                        {
                            if (k == AdressMin[1])
                            {
                                continue;
                            }
                            if (PriviousStepMatrix[AdressMin[0], k] != (double)Status.NotInit)
                                continue;
                            PriviousStepMatrix[AdressMin[0], k] = (double)Status.Empty;
                        }
                        i++;
                        step++;
                    }
                    CanGiveArray[AdressMin[0]] -= PriviousStepMatrix[AdressMin[0], AdressMin[1]];
                    NeedArray[AdressMin[1]] -= PriviousStepMatrix[AdressMin[0], AdressMin[1]];
                    if (step == inputArray.GetLength(0) + inputArray.GetLength(1) - 1) break;
                }
                if (step == inputArray.GetLength(0) + inputArray.GetLength(1) - 1) break;
            }
            //Console.WriteLine("Матрица поставок");
            //Console.Write(PriviousStepMatrix.ToString());
            return PriviousStepMatrix;
        }

        static public Matrix<double> BiuldFirstStep3(double[,] inputArray, double[] CanGiveArray/*N*/, double[] NeedArray/*M*/)
        {
            Matrix<double> PriviousStepMatrix = M.Dense(inputArray.GetLength(0), inputArray.GetLength(1), (double)Status.NotInit);
            int target = inputArray.GetLength(0) + inputArray.GetLength(1) - 1 - 2;
            Matrix<double> PriceTempMatrix = M.DenseOfArray(inputArray);
            int step = 0;
            while (true)
            {
                double[] diffColumn = new double[PriceTempMatrix.ColumnCount];
                for (int i = 0; i < PriceTempMatrix.ColumnCount; i++)
                {
                    if (PriceTempMatrix.Column(i).OrderBy(x => x).Where((x) =>
                    {
                        return x > 0;
                    }).ToArray().Count() >= 1)
                    {
                        double min = PriceTempMatrix.Column(i).OrderBy(x => x).ToArray()[0];
                        double nMin = PriceTempMatrix.Column(i).OrderBy(x => x).Where((x) => { return x > 0; }).ToArray()[0];
                        diffColumn[i] = nMin - min;
                    }
                }

                double[] diffRow = new double[PriceTempMatrix.RowCount];
                for (int i = 0; i < PriceTempMatrix.RowCount; i++)
                {
                    if (PriceTempMatrix.Row(i).OrderBy(x => x).Where((x) => { return x > 0; }).ToArray().Count() >= 1)
                    { 
                    

                    
                    double min = PriceTempMatrix.Row(i).OrderBy(x => x).ToArray()[0];
                    double nMin = PriceTempMatrix.Row(i).OrderBy(x => x).Where((x) => { return x > 0; }).ToArray()[0];
                    diffRow[i] = nMin - min;
                    }
                }

                int MinI = 0;
                int MinJ = 0;

                if (diffColumn.Max() > diffRow.Max())
                {
                    for (int i = 0; i < diffColumn.Length; i++)
                    {
                        if (diffColumn[i] == diffColumn.Max())
                        {
                            MinJ = i;
                        }
                    }
                    for (int i = 0; i < PriceTempMatrix.Column(MinJ).Count(); i++)
                    {
                        if (PriceTempMatrix.Column(MinJ)[i] == PriceTempMatrix.Column(MinJ).Where((x) => { return x > 0; }).Min())
                        {
                            MinI = i;
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < diffRow.Length; i++)
                    {
                        if (diffRow[i] == diffRow.Max())
                        {
                            MinI = i;
                        }
                    }
                    for (int i = 0; i < PriceTempMatrix.Row(MinI).Count(); i++)
                    {
                        if (PriceTempMatrix.Row(MinI)[i] == PriceTempMatrix.Row(MinI).Where((x) => { return x > 0; }).Min())
                        {
                            MinJ = i;
                        }
                    }
                }

                int[] AdressMin = new int[] { MinI, MinJ };
                /*if (step == target)
                {
                    AdressMin = new int[] { Max(CanGiveArray), Max(CanGiveArray) };
                }*/
                PriceTempMatrix[AdressMin[0], AdressMin[1]] = (double)Status.Empty;
                if (CanGiveArray[AdressMin[0]] >= NeedArray[AdressMin[1]])
                {
                    PriviousStepMatrix[AdressMin[0], AdressMin[1]] = NeedArray[AdressMin[1]];
                    for (int k = 0; k < CanGiveArray.Length; k++)
                    {
                        if (k == AdressMin[0])
                        {
                            continue;
                        }
                        if (PriviousStepMatrix[k, AdressMin[1]] != (double)Status.NotInit)
                            continue;
                        PriviousStepMatrix[k, AdressMin[1]] = (double)Status.Empty;
                        PriceTempMatrix[k, AdressMin[1]] = (double)Status.Empty;
                    }
                }
                else
                {
                    PriviousStepMatrix[AdressMin[0], AdressMin[1]] = CanGiveArray[AdressMin[0]];
                    for (int k = 0; k < NeedArray.Length; k++)
                    {
                        if (k == AdressMin[1])
                        {
                            continue;
                        }
                        if (PriviousStepMatrix[AdressMin[0], k] != (double)Status.NotInit)
                            continue;
                        PriviousStepMatrix[AdressMin[0], k] = (double)Status.Empty;
                        PriceTempMatrix[AdressMin[0], k] = (double)Status.Empty;
                    }
                }
                CanGiveArray[AdressMin[0]] -= PriviousStepMatrix[AdressMin[0], AdressMin[1]];
                NeedArray[AdressMin[1]] -= PriviousStepMatrix[AdressMin[0], AdressMin[1]];
                step++;
                if ((CanGiveArray.Sum() + NeedArray.Sum()) == 0)
                {
                    break;
                }
                
                //Console.Write(PriviousStepMatrix.ToString());
                //Console.Write(PriceTempMatrix.ToString());
            }
            return PriviousStepMatrix;
        }
    }
}
    
