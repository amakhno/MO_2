using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PotentialMethod
{
    class Program
    {
        static void Main(string[] args)
        {
            string[] input = File.ReadAllLines("input.txt");
            int M = input[0].Split(new char[] { ' ' }).Count()-1;
            int N = input.Length - 1;     //N x M
            double[,] inputArray = new double[N, M];
            for(int i = 0; i < N; i++)
            {
                string[] stolbec = input[i].Split(new char[] { ' ' });
                for (int j = 0; j < M; j++)
                {
                    inputArray[i, j] = Convert.ToDouble(stolbec[j]);
                }
            }
            double[] CanGiveArray = new double[N];
            double[] NeedArray = new double[M];


            for(int i = 0; i< N; i++)
            {
                CanGiveArray[i] = Convert.ToDouble(input[i].Split(new char[] { ' ' })[M]);
            }

            for (int i = 0; i < M; i++)
            {
                NeedArray[i] = Convert.ToDouble(input[N].Split(new char[] { ' ' })[i]);
            }
            if(CanGiveArray.Sum() != NeedArray.Sum())
            {
                Console.WriteLine("Плохие данные");
                Console.ReadKey();
                return;
            }

            var PriviousStepMatrix = PotencialMethod.BiuldFirstStep(inputArray, CanGiveArray, NeedArray);
            Console.WriteLine("Метод наименьшего");
            WriteAnswer(PriviousStepMatrix, input, N, M);
            for (int i = 0; i < N; i++)
            {
                CanGiveArray[i] = Convert.ToDouble(input[i].Split(new char[] { ' ' })[M]);
            }

            for (int i = 0; i < M; i++)
            {
                NeedArray[i] = Convert.ToDouble(input[N].Split(new char[] { ' ' })[i]);
            }
            PriviousStepMatrix = PotencialMethod.BiuldFirstStep2(inputArray, CanGiveArray, NeedArray);
            Console.WriteLine("\nМетод С-З угла");
            WriteAnswer(PriviousStepMatrix, input, N, M);
            for (int i = 0; i < N; i++)
            {
                CanGiveArray[i] = Convert.ToDouble(input[i].Split(new char[] { ' ' })[M]);
            }

            for (int i = 0; i < M; i++)
            {
                NeedArray[i] = Convert.ToDouble(input[N].Split(new char[] { ' ' })[i]);
            }
            PriviousStepMatrix = PotencialMethod.BiuldFirstStep3(inputArray, CanGiveArray, NeedArray);
            Console.WriteLine("\nМетод Фогеля");
            WriteAnswer(PriviousStepMatrix, input, N, M);
        }

        static public void WriteAnswer(Matrix<double> matrix, string[] input, int N, int M)
        {
            double[,] inputArray = new double[N, M];
            for (int i = 0; i < N; i++)
            {
                string[] stolbec = input[i].Split(new char[] { ' ' });
                for (int j = 0; j < M; j++)
                {
                    inputArray[i, j] = Convert.ToDouble(stolbec[j]);
                }
            }
            double[] CanGiveArray = new double[N];
            double[] NeedArray = new double[M];
            for (int i = 0; i < N; i++)
            {
                CanGiveArray[i] = Convert.ToDouble(input[i].Split(new char[] { ' ' })[M]);
            }

            for (int i = 0; i < M; i++)
            {
                NeedArray[i] = Convert.ToDouble(input[N].Split(new char[] { ' ' })[i]);
            }


            double result = 0;
            for(int i = 0; i<matrix.RowCount; i++)
            {
                for (int j = 0; j < matrix.ColumnCount; j++)
                {
                    if(matrix[i,j]!=-1)
                    {
                        Console.Write(String.Format("x[{0},{1}] = {2}, ", i, j, matrix[i, j]));
                        result += inputArray[i, j] * matrix[i, j];
                    }
                }                
            }
            Console.Write("Sum = " + result.ToString());
        } 
    }
}
