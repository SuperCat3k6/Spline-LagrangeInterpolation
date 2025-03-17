using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ScottPlot;

namespace Spline_LagrangeInterpolation
{
    class Program
    {
        static double Func(double x) => x - Math.Sin(x) - 0.25;

        static void Main()
        {
            int[] nodeCounts = { 5, 10, 20, 40 };
            int k = 1000; // Количество проверочных точек
            double xmin = 0, xmax = 2 * Math.PI;

            foreach (int n in nodeCounts)
            {
                double[] x = new double[n];
                double[] y = new double[n];

                // Заполняем узлы интерполяции
                for (int i = 0; i < n; i++)
                {
                    x[i] = xmin + i * (xmax - xmin) / (n - 1);
                    y[i] = Func(x[i]);
                }

                // Проверочные точки
                double[] xs = new double[k];
                double[] ysFunc = new double[k];
                double[] ysSpline = new double[k];
                double[] ysLagrange = new double[k];
                double[] errSpline = new double[k];
                double[] errLagrange = new double[k];

                for (int i = 0; i < k; i++)
                {
                    double xi = xmin + i * (xmax - xmin) / (k - 1);
                    ysFunc[i] = Func(xi);
                    ysSpline[i] = CubicSpline(x, y, xi);
                    ysLagrange[i] = LagrangeInterpolation(x, y, xi);
                    errSpline[i] = Math.Abs(ysFunc[i] - ysSpline[i]);
                    errLagrange[i] = Math.Abs(ysFunc[i] - ysLagrange[i]);
                    xs[i] = xi;
                }

                // Построение графиков погрешности
                PlotError(xs, errSpline, errLagrange, n);
            }
        }

        // Полином Лагранжа
        static double LagrangeInterpolation(double[] x, double[] y, double xVal)
        {
            double result = 0;
            for (int i = 0; i < x.Length; i++)
            {
                double term = y[i];
                for (int j = 0; j < x.Length; j++)
                {
                    if (j != i)
                        term *= (xVal - x[j]) / (x[i] - x[j]);
                }
                result += term;
            }
            return result;
        }

        // Кубический сплайн (натуральный)
        static double CubicSpline(double[] x, double[] y, double xVal)
        {
            int n = x.Length;
            double[] h = new double[n - 1];
            double[] alpha = new double[n - 1];
            double[] l = new double[n];
            double[] mu = new double[n];
            double[] z = new double[n];
            double[] c = new double[n];
            double[] b = new double[n - 1];
            double[] d = new double[n - 1];

            for (int i = 0; i < n - 1; i++)
                h[i] = x[i + 1] - x[i];

            for (int i = 1; i < n - 1; i++)
                alpha[i] = (3 / h[i]) * (y[i + 1] - y[i]) - (3 / h[i - 1]) * (y[i] - y[i - 1]);

            l[0] = 1;
            mu[0] = 0;
            z[0] = 0;

            for (int i = 1; i < n - 1; i++)
            {
                l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
                mu[i] = h[i] / l[i];
                z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
            }

            l[n - 1] = 1;
            z[n - 1] = 0;
            c[n - 1] = 0;

            for (int j = n - 2; j >= 0; j--)
            {
                c[j] = z[j] - mu[j] * c[j + 1];
                b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
                d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
            }

            int idx = FindSegment(x, xVal);
            double dx = xVal - x[idx];
            return y[idx] + b[idx] * dx + c[idx] * dx * dx + d[idx] * dx * dx * dx;
        }

        // Поиск сегмента
        static int FindSegment(double[] x, double xVal)
        {
            int i = 0, j = x.Length - 1;
            while (i < j)
            {
                int mid = (i + j) / 2;
                if (xVal > x[mid])
                    i = mid + 1;
                else
                    j = mid;
            }
            return Math.Max(0, i - 1);
        }

        // График ошибки
        static void PlotError(double[] xs, double[] errSpline, double[] errLagrange, int n)
        {

            var plt = new ScottPlot.Plot();
            var scatterSpline = plt.Add.Scatter(xs, errSpline);
            scatterSpline.Color = ScottPlot.Colors.Red;
            var scatterLagrange = plt.Add.Scatter(xs, errLagrange);
            scatterLagrange.Color = ScottPlot.Colors.Blue;
            scatterSpline.LegendText = "Ошибка кубического сплайна";
            scatterLagrange.LegendText = "Ошибка Лагранжа";
            plt.Title($"Распределение ошибки (n={n})");
            plt.Legend.Alignment = ScottPlot.Alignment.UpperLeft;
            plt.ShowLegend();

            string fileName = $"error_n{n}.png";
            plt.SavePng(fileName, 600, 400);
            Console.WriteLine($"График ошибки сохранен: {fileName}");
        }
    }
}