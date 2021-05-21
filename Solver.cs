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
            
            // ����� ���� �������� ������ � ������� ��������� ���������� ����!!!
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
             * ��� ������ �� ����������� �������, ������ � �������� ������� ������.
             * ����� �������, �� ������ �������� ��� ���������� ���������� ������� ��
             * �� �������� ������������� ������������� �������� � �� ����� ��� PreparedColumnsCount-� ������.
             *
             * ���� ��������������� �������� �� ��������, ������������� ��������� �� �������
             */

            // ������� ������.
            if (_preparedColumnsCount == Width) return;

            // ������: �� ���� �������� � ������� �� ������� (������� �������)
            // ����������� ���������. ������� �����
            if (rowIndex >= Height)
            {
                _preparedColumnsCount++;
                return;
            }

            // ����� m[r][c] == 1, ����� ������ ��� ������ �� m[r][c] 
            var divider = _matrix[rowIndex][columnIndex];
            if (divider == 0)
            {
                // �� ���� ������ ������, ���� ������� ��������� ��� ��������
                // (��������� ����� � �������)
                PrepareColumn(rowIndex + 1, columnIndex);
                return;
            }

            // ���� �������� ����������, �� �����.
            // m[r][*] / m[r][c], ����� m[r][c] == 1
            MultiplyLine(rowIndex, 1 / divider);

            // � �������� � � ������ ���-��� �� ���������.
            // m[..][c] - m[r] * m[..][c]
            for (var row = 0; row < Height; row++)
            {
                if (row == rowIndex) continue;
                AddMultipliedLine(row, rowIndex, -_matrix[row][columnIndex]);
            }

            // ������ m[*][c] == 0; m[r][c] == 1.
            // ������� ������� ������.

            // TODO TODO
            if (rowIndex != _preparedColumnsCount)
                SwitchLines(rowIndex, _preparedColumnsCount);
            _preparedColumnsCount++;
            _dependentVars[rowIndex, columnIndex] = true;
            _dependentVarsCount++;
        }

        public double[] Solve()
        {
            // �������� � ������������ ����
            // TODO: ������� ������ � ����������� ��������� ����������. ��� ���� ���������� �����
            // for (var dividerSearchStartIndex = 0; dividerSearchStartIndex < Height; dividerSearchStartIndex++) // ���� ������
            /*for (var dividerSearchStartIndex = 0; _preparedColumnsCount < Width; dividerSearchStartIndex++)
            {
                // �� �� ����� �������� ������ �������� ����� ��� �������������� �����.
                // �.�. ������ ������� �������������� ������ ���������� �����,
                // �� ��� ����� ����� �������� ������ ������ �� 1 ���� �� ������ ��������.
                
                // if (dividerSearchStartIndex == )

                PrepareColumn(dividerSearchStartIndex, _preparedColumnsCount);
            }*/


            while (_preparedColumnsCount < Width)
            {
                PrepareColumn(_dependentVarsCount, _preparedColumnsCount);
            }

            // ������ ������� ����� ����������� ���! ����� ��� ������� � ��������� �����, ������.

            // ������� ������������ ����, ��������� �� ���.
            // ���� ��� -- ����� ���������� ������ ��������� �������.
            for (int row = 0; row < Height; row++)
            {
                if (RowOnlyContainsZeros(row) && _freeMembers[row] != 0)
                    throw new NoSolutionException();
            }
            // ������ ������������� ���������� �������.
            
            var solution = new double[Width];
            var definedVars = new bool[Width];
            
            for (var row = Height - 1; row >= 0; row--)
            {
                if (RowOnlyContainsZeros(row)) continue;
                var encounteredDependentVar = false;
                for (var column = Width - 1; column >= 0 && !encounteredDependentVar; column--)
                {
                    // x_column -- ��������� ����������
                    // !!! TODO
                    // if (row == column && _dependentVars[column])
                    if (_dependentVars[row, column])
                    {
                        // ����������� ��� ��������� ����� �������� ����������� ����������
                        var sum = 0.0;
                        for (var i = column + 1; i < Width; i++)
                            sum += _matrix[row][i] * solution[i];
                        var tmp = (_freeMembers[row] - sum) / _matrix[row][column];
                        solution[column] = tmp;
                        
                        // ������� � ������ ������������, ����� �� ������� �� ��������� ����������
                        encounteredDependentVar = true;
                    }
                    else // x_column -- ��������� ����������
                    {
                        // ���� ���������� ���������� � �������, �� ������ �� ������.
                        // ���� �� �� ����������, �� ����������� �� ����� ��������
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