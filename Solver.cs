using System;
using System.Linq;

namespace GaussAlgorithm
{
    public class Solver
    {
        public double[] Solve(double[][] matrix, double[] freeMembers)
        {
            var system = new LinearEquationSystem(matrix, freeMembers);
            return system.Solve();
        }
    }

    public class LinearEquationSystem
    {
        private const double Epsilon = 1e-6;
        private readonly double[][] matrix;
        private readonly double[] freeMembers;
        public int Height => matrix.Length;
        public int Width => Height > 0 ? matrix[0].Length : 0;
        private int preparedColumnsCount;
        private readonly bool[][] dependentVars;
        private int dependentVarsCount;

        public LinearEquationSystem(double[][] matrix, double[] freeMembers)
        {
            this.matrix = matrix;
            this.freeMembers = freeMembers;
            
            dependentVars = new bool[Height][];
            for (var row = 0; row < Height; row++)
                dependentVars[row] = new bool [Width];
        }

        public void AddMultipliedLine(int resIndex, int addIndex, double multiplier)
        {
            for (var i = 0; i < Width; i++)
            {
                matrix[resIndex][i] += multiplier * matrix[addIndex][i];
                if (Math.Abs(matrix[resIndex][i]) < Epsilon) matrix[resIndex][i] = 0;
            }
            freeMembers[resIndex] += multiplier * freeMembers[addIndex];
            if (Math.Abs(freeMembers[resIndex]) < Epsilon) freeMembers[resIndex] = 0;
        }

        public void MultiplyLine(int resIndex, double multiplier)
        {
            matrix[resIndex] = matrix[resIndex].Select(x => x * multiplier).ToArray();
            freeMembers[resIndex] *= multiplier;
        }

        public void SwitchLines(int i, int j)
        {
            if (i < 0 || i >= Height || j < 0 || j >= Height) return;
            
            var tmpLine = matrix[i];
            matrix[i] = matrix[j];
            matrix[j] = tmpLine;

            var tmp = freeMembers[i];
            freeMembers[i] = freeMembers[j];
            freeMembers[j] = tmp;
            
            var tmpLine2 = dependentVars[i];
            dependentVars[i] = dependentVars[j];
            dependentVars[j] = tmpLine2;
        }

        private bool RowOnlyContainsZeros(int row) => matrix[row].All(x => x == 0);

        private void PrepareColumn(int rowIndex, int columnIndex)
        {
            /*
             * m[r][c] == 1, m[r][!= c] == 0
             * Как только мы подготовили столбец, строку с единицей двигаем наверх.
             * Таким образом, на каждой итерации при завершении подготовки столбца мы
             * мы начинаем рассматривать потенциальные делители с не менее чем PreparedColumnsCount-й строки.
             *
             * Если рассматриваемый делитель не подходит, рассматриваем следующий по столбцу
             */

            // Матрица готова.
            if (preparedColumnsCount == Width) return;

            // Случай: ни один делитель в столбце не подошёл (нулевой столбик)
            // Заканчиваем обработку. Столбик готов
            if (rowIndex >= Height)
            {
                preparedColumnsCount++;
                return;
            }

            // Чтобы m[r][c] == 1, будем делить всю строку на m[r][c] 
            var divider = matrix[rowIndex][columnIndex];
            if (divider == 0)
            {
                // На ноль делить нельзя, ищём другого кандидата для делителя
                // (следующее число в столбце)
                PrepareColumn(rowIndex + 1, columnIndex);
                return;
            }

            // Если делитель подходящий, то делим.
            // m[r][*] / m[r][c], чтобы m[r][c] == 1
            MultiplyLine(rowIndex, 1 / divider);

            // и вычитаем её в нужных кол-вах из остальных.
            // m[..][c] - m[r] * m[..][c]
            for (var row = 0; row < Height; row++)
            {
                if (row == rowIndex) continue;
                AddMultipliedLine(row, rowIndex, -matrix[row][columnIndex]);
            }

            // Теперь m[*][c] == 0; m[r][c] == 1.
            // Двигаем строчку наверх.
            
            dependentVars[rowIndex][columnIndex] = true;
            if (rowIndex != preparedColumnsCount)
                SwitchLines(rowIndex, preparedColumnsCount);
            preparedColumnsCount++;
            dependentVarsCount++;

        }

        public double[] Solve()
        {
            // Приводим к ступенчатому виду
            while (preparedColumnsCount < Width)
                PrepareColumn(dependentVarsCount, preparedColumnsCount);

            // Теперь матрица имеет ступенчатый вид! Можно это вынести в отдельный метод, кстати.

            // Провёдем исследование того, совместна ли она.
            // Если нет -- сразу возвращаем пустое множество решений.
            for (int row = 0; row < Height; row++)
            {
                if (RowOnlyContainsZeros(row) && freeMembers[row] != 0)
                    throw new NoSolutionException();
            }
            // Теперь рассматриваем совместную матрицу.
            
            var solution = new double[Width];
            var definedVars = new bool[Width];
            
            for (var row = Height - 1; row >= 0; row--)
            {
                if (RowOnlyContainsZeros(row)) continue;
                var encounteredDependentVar = false;
                for (var column = Width - 1; column >= 0 && !encounteredDependentVar; column--)
                {
                    // x_column -- зависимая переменная
                    // !!! TODO
                    // if (row == column && _dependentVars[column])
                    if (dependentVars[row][column])
                    {
                        // подставляем все найденные ранее значения независимых переменных
                        var sum = 0.0;
                        for (var i = column + 1; i < Width; i++)
                            sum += matrix[row][i] * solution[i];
                        var tmp = (freeMembers[row] - sum) / matrix[row][column];
                        solution[column] = tmp;
                        
                        // подсчёт в строке оканчивается, когда мы доходим до зависимой переменной
                        encounteredDependentVar = true;
                    }
                    else // x_column -- свободная переменная
                    {
                        // Если переменная определена в решении, то ничего не делаем.
                        // Если же не определена, то подставляем ей любое значение
                        if (!definedVars[column]) 
                            solution[column] = 0;
                    }
                    
                    definedVars[column] = true;
                }
            }

            return solution;
        }
    }
}