using System;
using System.Linq;

namespace GaussAlgorithm
{
    public class Solver
    {
        // TODO: Implement LINQ!!!

        public double[] Solve(double[][] matrix, double[] freeMembers)
        {
            var system = new LinearEquationSystem(matrix, freeMembers);
            return system.Solve();
        }
    }

    partial class NoSolutionException : Exception
    {
        public NoSolutionException() {}
    }

    public class LinearEquationSystem
    {
        private const double Epsilon = 1e-6;
        private readonly double[][] _matrix;
        private readonly double[] _freeMembers;
        public int Height => _matrix.Length;
        public int Width => Height > 0 ? _matrix[0].Length : 0;
        private int _preparedColumnsCount;
        // private readonly bool[] _dependentVars;
        private readonly bool[,] _dependentVars;
        // private int _dependentVarsCount => _dependentVars.Count(x => x);
        // private bool _isStepMatrix => _preparedColumnsCount == Width;
        private int _dependentVarsCount;

        public LinearEquationSystem(double[][] matrix, double[] freeMembers)
        {
            this._matrix = matrix;
            this._freeMembers = freeMembers;
            _dependentVars = new bool[Height, Width];
        }

        public void AddMultipliedLine(int resIndex, int addIndex, double multiplier)
        {
            for (var i = 0; i < Width; i++)
            {
                _matrix[resIndex][i] += multiplier * _matrix[addIndex][i];
                if (Math.Abs(_matrix[resIndex][i]) < Epsilon) _matrix[resIndex][i] = 0;
            }
            _freeMembers[resIndex] += multiplier * _freeMembers[addIndex];
            if (Math.Abs(_freeMembers[resIndex]) < Epsilon) _freeMembers[resIndex] = 0;

            Console.WriteLine("Add {0}th line multiplied by {1} to line {2}", addIndex + 1, multiplier, resIndex + 1);
            ViewMatrix();
        }

        public void MultiplyLine(int resIndex, double multiplier)
        {
            _matrix[resIndex] = _matrix[resIndex].Select(x => x * multiplier).ToArray();
            _freeMembers[resIndex] *= multiplier;

            Console.WriteLine("Multiply line {0} by {1}", resIndex + 1, multiplier);
            ViewMatrix();
        }

        public void SwitchLines(int i, int j)
        {
            if (i < 0 || i >= Height || j < 0 || j >= Height) return;
            
            var tmpLine = _matrix[i];
            _matrix[i] = _matrix[j];
            _matrix[j] = tmpLine;

            var tmp = _freeMembers[i];
            _freeMembers[i] = _freeMembers[j];
            _freeMembers[j] = tmp;
            
            // здесь надо поменять строки в матрице зависимых переменных тоже!!!
            /*var tmpLine2 = _matrix[i];
            _matrix[i] = _matrix[j];
            _matrix[j] = tmpLine2;*/

            Console.WriteLine("Switched {0} and {1} line", i, j);
            ViewMatrix();
        }

        public void ViewMatrix()
        {
            for (var row = 0; row < Height; row++)
            {
                foreach (var e in _matrix[row])
                    Console.Write(e + " ");
                Console.Write("| " + _freeMembers[row] + "\n");
            }
        }

        private bool RowOnlyContainsZeros(int row) => _matrix[row].All(x => x == 0);

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
            if (_preparedColumnsCount == Width) return;

            // Случай: ни один делитель в столбце не подошёл (нулевой столбик)
            // Заканчиваем обработку. Столбик готов
            if (rowIndex >= Height)
            {
                _preparedColumnsCount++;
                return;
            }

            // Чтобы m[r][c] == 1, будем делить всю строку на m[r][c] 
            var divider = _matrix[rowIndex][columnIndex];
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
                AddMultipliedLine(row, rowIndex, -_matrix[row][columnIndex]);
            }

            // Теперь m[*][c] == 0; m[r][c] == 1.
            // Двигаем строчку наверх.

            // TODO TODO
            if (rowIndex != _preparedColumnsCount)
                SwitchLines(rowIndex, _preparedColumnsCount);
            _preparedColumnsCount++;
            _dependentVars[rowIndex, columnIndex] = true;
            _dependentVarsCount++;
        }

        public double[] Solve()
        {
            // Приводим к ступенчатому виду
            // TODO: главная ошибка в определении зависимых переменных. Это надо исправлять ЗДЕСЬ
            // for (var dividerSearchStartIndex = 0; dividerSearchStartIndex < Height; dividerSearchStartIndex++) // было раньше
            /*for (var dividerSearchStartIndex = 0; _preparedColumnsCount < Width; dividerSearchStartIndex++)
            {
                // Мы не хотим начинать искать делитель среди уже подготовленных строк.
                // Т.к. каждая успешно подготовленная строка сдвигается вверх,
                // то для этого можно сдвигать начало поиска на 1 вниз на каждой итерации.
                
                // if (dividerSearchStartIndex == )

                PrepareColumn(dividerSearchStartIndex, _preparedColumnsCount);
            }*/


            while (_preparedColumnsCount < Width)
            {
                PrepareColumn(_dependentVarsCount, _preparedColumnsCount);
            }

            // Теперь матрица имеет ступенчатый вид! Можно это вынести в отдельный метод, кстати.

            // Провёдем исследование того, совместна ли она.
            // Если нет -- сразу возвращаем пустое множество решений.
            for (int row = 0; row < Height; row++)
            {
                if (RowOnlyContainsZeros(row) && _freeMembers[row] != 0)
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
                    if (_dependentVars[row, column])
                    {
                        // подставляем все найденные ранее значения независимых переменных
                        var sum = 0.0;
                        for (var i = column + 1; i < Width; i++)
                            sum += _matrix[row][i] * solution[i];
                        var tmp = (_freeMembers[row] - sum) / _matrix[row][column];
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