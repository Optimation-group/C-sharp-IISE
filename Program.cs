using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ILOG.Concert;
using ILOG.CPLEX;
using System.IO;
namespace 模型7._11
{
    class Program
    {
        public const int customer = 9;
        public const int Numtruck = 1;
        public const int rongl = 50;
        public long M = customer + 10;
        const int N = 10000;
        public const double speed_D = 48;
        public const double speed_T = 30;
        public const double ef = 0.5;
        public const double mk = 200;
        public const double md = 2.3;
        public const double ck = 0.693 * speed_T;
        public const double cd = 0.19 * speed_D;
        public const double cw = ck / 3;
        public const double cg = 100;
        public const double g = 15;

        public int r_n = customer;
        public int numtruck = Numtruck;

        public double[,] rloc_xy;
        public double[] q;
        public double[] efj;

        public double[,] timeD;
        public double[,,] timeT;

        void input()
        {

            rloc_xy = new double[r_n + 2, 2];
            timeT = new double[r_n + 2, r_n + 2, numtruck];
            timeD = new double[r_n + 2, r_n + 2];
            q = new double[r_n + 2];
            efj = new double[r_n + 2];
            string di;

            StreamReader WWW = new StreamReader(@"F:\Input\9-1.txt", false);

            for (int j = 0; j < r_n + 2; j++)
            {
                di = WWW.ReadLine();
                rloc_xy[j, 0] = Math.Round(Convert.ToSingle(di), 2);
                di = WWW.ReadLine();
                rloc_xy[j, 1] = Math.Round(Convert.ToSingle(di), 2);

            }
            for (int j = 0; j < r_n + 2; j++)
            {
                di = WWW.ReadLine();
                q[j] = Math.Round(Convert.ToSingle(di), 2);

            }
            for (int i = 1; i < r_n + 1; i++)
            {
                di = WWW.ReadLine();
                efj[i] = Math.Round(Convert.ToSingle(di), 2);

            }
            for (int i = 0; i < r_n + 2; i++)
            {
                for (int j = 0; j < r_n + 2; j++)
                {
                    for (int k = 0; k < numtruck; k++)
                    {
                        di = WWW.ReadLine();
                        timeT[i, j, k] = Math.Round(Convert.ToSingle(di), 2);

                    }
                }
            }
            for (int i = 0; i < r_n + 2; i++)
            {
                for (int j = 0; j < r_n + 2; j++)
                {
                    di = WWW.ReadLine();
                    timeD[i, j] = Math.Round(Convert.ToSingle(di), 2);

                }
            }

        }
        void solve(ref double OBBBBJ)
        {
            Cplex Model = new Cplex();
            INumVar[][][] A = new INumVar[r_n + 2][][];
            for (int i = 0; i < r_n + 2; i++)
            {
                A[i] = new INumVar[r_n + 2][];
                for (int j = 0; j < r_n + 2; j++)
                {
                    A[i][j] = new INumVar[numtruck];
                    for (int k = 0; k < numtruck; k++)
                    {
                        if (i == j)
                            A[i][j][k] = Model.NumVar(0, 0, NumVarType.Bool, "alfa_i" + i.ToString() + "_j" + j.ToString() + "_k" + k.ToString());
                        if (i != j)
                            A[i][j][k] = Model.NumVar(0, 1, NumVarType.Bool, "alfa_i" + i.ToString() + "_j" + j.ToString() + "_k" + k.ToString());
                    }
                }
            }
            INumVar[][] a_T = new INumVar[r_n + 2][];
            for (int i = 0; i < r_n + 2; i++)
            {
                a_T[i] = new INumVar[numtruck];
                for (int k = 0; k < numtruck; k++)
                {
                    a_T[i][k] = Model.NumVar(0, double.MaxValue, NumVarType.Float, "a_T_i" + i.ToString() + "_k" + k.ToString());
                }
            }
            INumVar[][][][] beta_ijg = new INumVar[r_n + 2][][][];
            for (int i = 0; i < r_n + 2; i++)
            {
                beta_ijg[i] = new INumVar[r_n + 2][][];
                for (int j = 0; j < r_n + 1; j++)
                {
                    beta_ijg[i][j] = new INumVar[r_n + 2][];
                    for (int h = 0; h < r_n + 2; h++)
                    {
                        beta_ijg[i][j][h] = new INumVar[numtruck];
                        for (int d = 0; d < numtruck; d++)
                        {
                            if (i != j && i != h && j != h)
                                beta_ijg[i][j][h][d] = Model.NumVar(0, 1, NumVarType.Bool, "beta_i" + i.ToString() + "_j" + j.ToString() + "_h" + h.ToString() + "_d" + d.ToString());
                            else
                                beta_ijg[i][j][h][d] = Model.NumVar(0, 0, NumVarType.Bool, "beta_i" + i.ToString() + "_j" + j.ToString() + "_h" + h.ToString() + "_d" + d.ToString());
                        }
                    }
                }
            }

            INumVar[][] a_D = new INumVar[r_n + 2][];
            for (int i = 0; i < r_n + 2; i++)
            {
                a_D[i] = new INumVar[numtruck];
                for (int k = 0; k < numtruck; k++)
                {
                    a_D[i][k] = Model.NumVar(0, double.MaxValue, NumVarType.Float, "a_D_i" + i.ToString() + "_k" + k.ToString());
                }
            }

            INumVar[][] U = new INumVar[r_n + 2][];
            for (int i = 0; i < r_n + 2; i++)
            {
                U[i] = new INumVar[numtruck];
                for (int k = 0; k < numtruck; k++)
                {
                    U[i][k] = Model.NumVar(0, r_n + 2, NumVarType.Int, "U_i" + i.ToString() + "_k" + k.ToString());
                }
            }
            INumVar[] E = new INumVar[numtruck];
            for (int k = 0; k < numtruck; k++)
            {
                E[k] = Model.NumVar(0, 1, NumVarType.Bool, "E_k" + k.ToString());
            }

            INumVar[][][] p_ijk = new INumVar[r_n + 2][][];
            for (int i = 0; i < r_n + 2; i++)
            {
                p_ijk[i] = new INumVar[r_n + 2][];
                for (int j = 0; j < r_n + 2; j++)
                {
                    p_ijk[i][j] = new INumVar[numtruck];
                    for (int k = 0; k < numtruck; k++)
                    {
                        p_ijk[i][j][k] = Model.NumVar(0, 1, NumVarType.Bool, "p_i" + i.ToString() + "_j" + j.ToString() + "_k" + k.ToString());
                    }
                }
            }

            INumVar[][] X = new INumVar[r_n + 2][];
            for (int i = 0; i < r_n + 2; i++)
            {
                X[i] = new INumVar[numtruck];
                for (int k = 0; k < numtruck; k++)
                {
                    X[i][k] = Model.NumVar(0, double.MaxValue, NumVarType.Float, "X_i" + i.ToString() + "_k" + k.ToString());
                }
            }

            #region Constraint 
            INumExpr[] expr1 = new INumExpr[1];
            IRange[] rang1 = new IRange[numtruck];
            for (int k = 0; k < numtruck; k++)
            {
                expr1[0] = A[0][0][0];
                expr1[0] = Model.Sum(expr1[0], Model.Prod(-1, A[0][0][0]));
                for (int j = 1; j < r_n + 2; j++)
                {
                    expr1[0] = Model.Sum(expr1[0], A[0][j][k]);
                }
                rang1[k] = Model.AddRange(1, 1);
                rang1[k].Expr = expr1[0];
            }
            #endregion

            #region Constraint 
            INumExpr[] expr2 = new INumExpr[1];
            IRange[] rang2 = new IRange[numtruck];
            for (int k = 0; k < numtruck; k++)
            {
                expr2[0] = A[0][r_n + 1][0];
                expr2[0] = Model.Sum(expr2[0], Model.Prod(-1, A[0][r_n + 1][0]));
                for (int i = 0; i < r_n + 1; i++)
                {
                    expr2[0] = Model.Sum(expr2[0], A[i][r_n + 1][k]);
                }

                rang2[k] = Model.AddRange(1, 1);
                rang2[k].Expr = expr2[0];
            }

            #endregion

            #region Constraint 

            for (int j = 1; j < r_n + 1; j++)
            {
                for (int k = 0; k < numtruck; k++)
                {
                    INumExpr[] expr3_1 = new INumExpr[1];
                    INumExpr[] expr3_2 = new INumExpr[1];
                    expr3_1[0] = A[0][0][0];
                    expr3_1[0] = Model.Sum(expr3_1[0], Model.Prod(-1, A[0][0][0]));
                    expr3_2[0] = A[0][0][0];
                    expr3_2[0] = Model.Sum(expr3_2[0], Model.Prod(-1, A[0][0][0]));
                    for (int i = 0; i < r_n + 1; i++)
                    {
                        if (i != j)
                        {
                            expr3_1[0] = Model.Sum(expr3_1[0], A[i][j][k]);
                        }
                    }
                    for (int h = 1; h < r_n + 2; h++)
                    {
                        if (h != j)
                        {
                            expr3_2[0] = Model.Sum(expr3_2[0], A[j][h][k]);
                        }
                    }
                    Model.AddLe(expr3_1[0], 1);
                    Model.AddLe(expr3_2[0], 1);
                    Model.AddEq(Model.Sum(expr3_2[0], Model.Prod(-1, expr3_1[0])), 0);
                }
            }
            #endregion

            #region Constraint
            INumExpr[] expr4 = new INumExpr[1];
            IRange[] rang4 = new IRange[r_n + 1];
            for (int j = 1; j < r_n + 1; j++)
            {
                expr4[0] = beta_ijg[0][1][2][0];
                expr4[0] = Model.Sum(expr4[0], Model.Prod(-1, beta_ijg[0][1][2][0]));
                for (int d = 0; d < numtruck; d++)
                {
                    for (int i = 0; i < r_n + 1; i++)
                    {
                        for (int h = 1; h < r_n + 2; h++)
                        {
                            if (i != j && j != h & i != h)
                            {
                                expr4[0] = Model.Sum(expr4[0], beta_ijg[i][j][h][d]);
                            }
                        }
                    }
                }

                rang4[j] = Model.AddRange(0, 1);
                rang4[j].Expr = expr4[0];
            }
            #endregion

            #region Constraint
            INumExpr[] expr5 = new INumExpr[1];
            IRange[] rang5 = new IRange[r_n + 1];
            for (int i = 0; i < r_n + 1; i++)
            {
                expr5[0] = beta_ijg[0][1][2][0];
                expr5[0] = Model.Sum(expr5[0], Model.Prod(-1, beta_ijg[0][1][2][0]));
                for (int d = 0; d < numtruck; d++)
                {
                    for (int h = 1; h < r_n + 2; h++)
                    {
                        for (int j = 1; j < r_n + 1; j++)
                        {
                            if (i != j && j != h && i != h)
                            {
                                expr5[0] = Model.Sum(expr5[0], beta_ijg[i][j][h][d]);
                            }
                        }
                    }
                }

                rang5[i] = Model.AddRange(0, 1);
                rang5[i].Expr = expr5[0];
            }
            #endregion

            #region Constraint
            INumExpr[] expr6 = new INumExpr[1];
            IRange[] rang6 = new IRange[r_n + 2];
            for (int h = 1; h < r_n + 2; h++)
            {
                for (int d = 0; d < numtruck; d++)
                {
                    expr6[0] = beta_ijg[0][1][2][0];
                    expr6[0] = Model.Sum(expr6[0], Model.Prod(-1, beta_ijg[0][1][2][0]));
                    for (int i = 0; i < r_n + 1; i++)
                    {

                        for (int j = 1; j < r_n + 1; j++)
                        {
                            if (i != j && j != h && i != h)
                            {
                                expr6[0] = Model.Sum(expr6[0], beta_ijg[i][j][h][d]);
                            }
                        }
                    }

                    rang6[h] = Model.AddRange(0, 1);
                    rang6[h].Expr = expr6[0];
                }

            }
            #endregion

            #region Constraint       
            INumExpr[] expr7_1 = new INumExpr[1];
            INumExpr[] expr7_2 = new INumExpr[1];
            IRange[] rang7 = new IRange[r_n + 1];
            for (int j = 1; j < r_n + 1; j++)
            {
                expr7_1[0] = A[0][1][0];
                expr7_1[0] = Model.Sum(expr7_1[0], Model.Prod(-1, A[0][1][0]));
                expr7_2[0] = beta_ijg[0][1][2][0];
                expr7_2[0] = Model.Sum(expr7_2[0], Model.Prod(-1, beta_ijg[0][1][2][0]));
                for (int k = 0; k < numtruck; k++)
                {
                    for (int i = 0; i < r_n + 1; i++)
                    {
                        if (i != j)
                        {
                            expr7_1[0] = Model.Sum(expr7_1[0], A[i][j][k]);
                        }
                    }
                }


                for (int d = 0; d < numtruck; d++)
                {
                    for (int h = 1; h < r_n + 2; h++)
                    {
                        for (int i = 0; i < r_n + 1; i++)
                        {
                            if (i != j && i != h && j != h)
                            {
                                expr7_2[0] = Model.Sum(expr7_2[0], beta_ijg[i][j][h][d]);
                            }
                        }
                    }
                }
                rang7[j] = Model.AddRange(1, 1);
                rang7[j].Expr = Model.Sum(expr7_1[0], expr7_2[0]);
            }
            #endregion

            #region Constraint 
            INumExpr[] expr8_1 = new INumExpr[1];
            INumExpr[] expr8_2 = new INumExpr[1];
            for (int i = 1; i < r_n + 1; i++)
            {
                for (int j = 1; j < r_n + 1; j++)
                {
                    for (int h = 1; h < r_n + 2; h++)
                    {
                        if (i != j && i != h && j != h)
                        {
                            for (int k = 0; k < numtruck; k++)
                            {
                                expr8_1[0] = A[0][1][0];
                                expr8_1[0] = Model.Sum(expr8_1[0], Model.Prod(-1, A[0][1][0]));
                                for (int j1 = 0; j1 < r_n + 1; j1++)
                                {
                                    if (j1 != i)
                                    {
                                        expr8_1[0] = Model.Sum(expr8_1[0], A[j1][i][k]);
                                    }
                                }
                                expr8_2[0] = A[0][1][0];
                                expr8_2[0] = Model.Sum(expr8_2[0], Model.Prod(-1, A[0][1][0]));
                                for (int i1 = 1; i1 < r_n + 1; i1++)
                                {
                                    if (i1 != h)
                                    {
                                        expr8_2[0] = Model.Sum(expr8_2[0], A[i1][h][k]);
                                    }
                                }
                                Model.AddLe(Model.Prod(2, beta_ijg[i][j][h][k]), Model.Sum(expr8_2[0], expr8_1[0]));
                            }
                        }
                    }
                }
            }
            #endregion

            #region Constraint 
            INumExpr[] expr9 = new INumExpr[1];
            for (int k = 0; k < numtruck; k++)
            {
                for (int j = 1; j < r_n + 1; j++)
                {
                    for (int h = 1; h < r_n + 2; h++)
                    {
                        if (j != h)
                        {
                            expr9[0] = A[0][1][0];
                            expr9[0] = Model.Sum(expr9[0], Model.Prod(-1, A[0][1][0]));
                            for (int i = 0; i < r_n + 1; i++)
                            {
                                if (i != h)
                                {
                                    expr9[0] = Model.Sum(expr9[0], A[i][h][k]);
                                }
                            }
                            Model.AddLe(beta_ijg[0][j][h][k], expr9[0]);
                        }
                    }
                }
            }
            #endregion

            #region Constraint    
            for (int k = 0; k < numtruck; k++)
            {
                for (int i = 0; i < r_n + 1; i++)
                {
                    for (int j = 1; j < r_n + 2; j++)
                    {
                        if (i != j)
                        {
                            Model.AddGe(Model.Sum(a_T[j][k], Model.Prod(M, Model.Sum(1, Model.Prod(-1, A[i][j][k])))), Model.Sum(X[i][k], timeT[i, j, k]));
                            //Model.AddLe(Model.Sum(a_T[i][k], timeT[i, j, k]), expr20_2[0]);
                        }
                    }
                }
            }
            #endregion

            #region Constraint 
            INumExpr[] expr21 = new INumExpr[1];
            for (int j = 1; j < r_n + 1; j++)
            {
                for (int i = 0; i < r_n + 1; i++)
                {
                    if (i != j)
                    {
                        for (int d = 0; d < numtruck; d++)
                        {
                            expr21[0] = beta_ijg[0][0][0][0];
                            expr21[0] = Model.Sum(expr21[0], Model.Prod(-1, beta_ijg[0][0][0][0]));
                            for (int h = 1; h < r_n + 2; h++)
                            {
                                if (i != h && j != h)
                                {
                                    expr21[0] = Model.Sum(expr21[0], beta_ijg[i][j][h][d]);
                                }
                            }
                            expr21[0] = Model.Prod(M, Model.Sum(1, Model.Prod(-1, expr21[0])));
                            Model.AddGe(Model.Sum(expr21[0], a_D[j][d]), Model.Sum(timeD[i, j], X[i][d]));
                        }
                    }
                }
            }
            #endregion

            #region Constraint
            INumExpr[] expr22 = new INumExpr[2];

            for (int h = 1; h < r_n + 2; h++)
            {
                for (int j = 1; j < r_n + 1; j++)
                {
                    if (h != j)
                    {
                        for (int d = 0; d < numtruck; d++)
                        {
                            expr22[0] = beta_ijg[0][1][2][0];
                            expr22[0] = Model.Sum(expr22[0], Model.Prod(-1, beta_ijg[0][1][2][0]));
                            for (int i = 0; i < r_n + 1; i++)
                            {
                                if (i != h && j != i)
                                {
                                    expr22[0] = Model.Sum(expr22[0], beta_ijg[i][j][h][d]);
                                }
                            }
                            expr22[0] = Model.Prod(1 * M, Model.Sum(1, Model.Prod(-1, expr22[0])));
                            expr22[0] = Model.Sum(expr22[0], a_D[h][d], Model.Prod(-1, X[j][d]));
                            Model.AddGe(expr22[0], timeD[j, h]);
                        }
                    }
                }
            }
            #endregion

            #region Constraint              
            INumExpr[] expr14 = new INumExpr[1];
            for (int k = 0; k < numtruck; k++)
            {
                for (int i = 1; i < r_n + 1; i++)
                {
                    for (int h = 1; h < r_n + 2; h++)
                    {
                        if (h != i)
                        {
                            expr14[0] = beta_ijg[1][2][3][0];
                            expr14[0] = Model.Sum(expr14[0], Model.Prod(-1, beta_ijg[1][2][3][0]));
                            for (int j = 1; j < r_n + 1; j++)
                            {
                                if (i != j && j != h)
                                {
                                    expr14[0] = Model.Sum(expr14[0], beta_ijg[i][j][h][k]);
                                }
                            }
                            expr14[0] = Model.Prod(r_n + 2, Model.Sum(1, Model.Prod(-1, expr14[0])));
                            Model.AddLe(Model.Sum(1, U[i][k]), Model.Sum(expr14[0], U[h][k]));
                        }
                    }
                }
            }
            #endregion

            #region Constraint [15]       
            for (int k = 0; k < numtruck; k++)
            {
                for (int i = 0; i < r_n + 1; i++)
                {
                    for (int j = 1; j < r_n + 2; j++)
                    {
                        if (i != j)
                        {
                            Model.AddLe(Model.Sum(U[i][k], 1), Model.Sum(U[j][k], Model.Prod(r_n + 2, Model.Sum(1, Model.Prod(-1, A[i][j][k])))));
                            //Model.AddGe(Model.Sum(U[i][k], 1), Model.Sum(U[j][k], Model.Prod((r_n + 2), Model.Sum(1, Model.Prod(-1, A[i][j][k])))));
                        }
                    }
                }
            }
            #endregion

            #region Constraint 
            INumExpr[] expr17 = new INumExpr[1];
            for (int k = 0; k < numtruck; k++)
            {
                for (int j = 1; j < r_n + 2; j++)
                {
                    expr17[0] = beta_ijg[0][1][2][0];
                    expr17[0] = Model.Sum(expr17[0], Model.Prod(beta_ijg[0][1][2][0], -1));
                    for (int i = 0; i < r_n + 1; i++)
                    {
                        if (i != j)
                            expr17[0] = Model.Sum(expr17[0], Model.Prod(r_n + 2, A[i][j][k]));
                    }
                    Model.AddLe(U[j][k], expr17[0]);
                }

            }
            #endregion

            #region Constraint 
            INumExpr[] expr17_1 = new INumExpr[1];

            for (int k = 0; k < numtruck; k++)
            {
                for (int i = 1; i < r_n + 1; i++)
                {
                    for (int j = 1; j < r_n + 1; j++)
                    {
                        if (i != j)
                        {
                            Model.AddLe(U[j][k], Model.Sum(U[i][k], Model.Prod(r_n + 2, p_ijk[i][j][k])));
                            //Model.AddLe(Model.Sum(U[i][k], 1), Model.Sum(U[j][k], Model.Prod((r_n + 2), Model.Sum(1, Model.Prod(-1, p_ijk[i][j][k])))));
                        }
                    }
                }
            }
            #endregion

            #region Constraint

            for (int k = 0; k < numtruck; k++)
            {
                for (int i = 1; i < r_n + 1; i++)
                {
                    for (int j = 1; j < r_n + 1; j++)
                    {
                        if (i != j)
                        {
                            Model.AddLe(Model.Sum(U[i][k], 1), Model.Sum(U[j][k], Model.Prod(r_n + 2, Model.Sum(1, Model.Prod(-1, p_ijk[i][j][k])))));
                        }
                    }
                }
            }
            #endregion


            #region Constraint           
            for (int i1 = 0; i1 < r_n + 1; i1++)
            {
                for (int h1 = 1; h1 < r_n + 2; h1++)
                {
                    for (int i2 = 1; i2 < r_n + 1; i2++)
                    {
                        for (int k = 0; k < numtruck; k++)
                        {
                            if (h1 != i1 && i2 != i1 && i2 != h1)
                            {
                                INumExpr[] expr19_1 = new INumExpr[1];
                                INumExpr[] expr19_2 = new INumExpr[1];
                                expr19_1[0] = beta_ijg[0][1][2][0];
                                expr19_1[0] = Model.Sum(expr19_1[0], Model.Prod(beta_ijg[0][1][2][0], -1));
                                for (int j1 = 1; j1 < r_n + 1; j1++)
                                {
                                    if (j1 != i2 && i1 != j1 && j1 != h1)
                                    {
                                        expr19_1[0] = Model.Sum(expr19_1[0], beta_ijg[i1][j1][h1][k]);
                                    }
                                }
                                expr19_2[0] = beta_ijg[0][1][2][0];
                                expr19_2[0] = Model.Sum(expr19_2[0], Model.Prod(beta_ijg[0][1][2][0], -1));
                                for (int j2 = 1; j2 < r_n + 1; j2++)
                                {
                                    for (int h2 = 1; h2 < r_n + 2; h2++)
                                    {
                                        if (j2 != h1 && j2 != i1 && j2 != i2 && h2 != i1 && h2 != h1 && i2 != h2 && j2 != h2)
                                        {
                                            expr19_2[0] = Model.Sum(expr19_2[0], beta_ijg[i2][j2][h2][k]);
                                        }
                                    }
                                }
                                Model.AddGe(X[i2][k], Model.Sum(a_D[h1][k], Model.Prod(-M, Model.Sum(3, Model.Prod(-1, Model.Sum(Model.Sum(expr19_1[0], expr19_2[0]), p_ijk[i1][i2][k]))))));

                            }
                        }
                    }
                }

            }
            #endregion

            #region Constraint          
            for (int i = 0; i < r_n + 1; i++)
            {
                for (int j = 1; j < r_n + 1; j++)
                {
                    for (int h = 1; h < r_n + 2; h++)
                    {
                        if (i != j && j != h && i != h)
                        {
                            for (int d = 0; d < numtruck; d++)
                            {
                                Model.AddLe(Model.Prod(timeD[i, j] + timeD[j, h], beta_ijg[i][j][h][d]), efj[j]);
                            }
                        }
                    }
                }
            }
            #endregion

            #region Constraint          
            for (int i = 0; i < r_n + 1; i++)
            {
                for (int j = 1; j < r_n + 1; j++)
                {
                    for (int h = 1; h < r_n + 2; h++)
                    {
                        if (i != j && j != h && i != h)
                        {
                            for (int d = 0; d < numtruck; d++)
                            {
                                Model.AddLe(Model.Sum(X[h][d], Model.Prod(-1, X[i][d])), Model.Sum(Model.Prod(M, Model.Sum(1, Model.Prod(-1, beta_ijg[i][j][h][d]))), efj[j]));
                                Model.AddLe(Model.Sum(Model.Sum(a_D[h][d], timeD[i, j]), Model.Prod(-1, X[j][d])), Model.Sum(Model.Prod(M, Model.Sum(1, Model.Prod(-1, beta_ijg[i][j][h][d]))), efj[j]));
                            }
                        }
                    }
                }
            }
            #endregion

            #region Constraint 
            for (int k = 0; k < numtruck; k++)
            {
                INumExpr[] expr24_1 = new INumExpr[1];
                INumExpr[] expr24_2 = new INumExpr[1];
                expr24_1[0] = A[0][1][0];
                expr24_1[0] = Model.Sum(expr24_1[0], Model.Prod(-1, A[0][1][0]));
                for (int i = 1; i < r_n + 1; i++)
                {
                    for (int j = 0; j < r_n + 1; j++)
                    {
                        expr24_1[0] = Model.Sum(expr24_1[0], Model.Prod(q[i], A[j][i][k]));
                    }
                }
                expr24_2[0] = Model.Prod(q[1], beta_ijg[0][1][2][0]);
                expr24_2[0] = Model.Sum(expr24_2[0], Model.Prod(-1, beta_ijg[0][1][2][0]));
                for (int i = 1; i < r_n + 1; i++)
                {
                    for (int j = 0; j < r_n + 1; j++)
                    {
                        for (int h = 1; h < r_n + 2; h++)
                        {
                            if (i != j && h != i && h != j)
                            {
                                expr24_2[0] = Model.Sum(expr24_2[0], Model.Prod(q[i], beta_ijg[j][i][h][k]));
                            }
                        }
                    }
                }
                Model.AddLe(Model.Sum(expr24_1[0], expr24_2[0]), mk);
            }
            #endregion

            #region Constraint   
            INumExpr[] expr25 = new INumExpr[1];
            for (int d = 0; d < numtruck; d++)
            {
                for (int i = 1; i < r_n + 1; i++)
                {
                    expr25[0] = beta_ijg[0][0][0][0];
                    expr25[0] = Model.Sum(expr25[0], Model.Prod(-1, beta_ijg[0][0][0][0]));
                    for (int j = 0; j < r_n + 1; j++)
                    {
                        for (int h = 1; h < r_n + 2; h++)
                        {
                            if (i != j && h != i && h != j)
                            {
                                expr25[0] = Model.Sum(expr25[0], beta_ijg[j][i][h][d]);
                            }
                        }
                    }
                    Model.AddLe(Model.Prod(q[i], expr25[0]), md);
                }
            }
            #endregion

            #region   Constraint 
            for (int i = 0; i < r_n + 1; i++)
            {
                for (int j = 1; j < r_n + 1; j++)
                {
                    for (int k = 0; k < numtruck; k++)
                    {
                        if (i != j)
                            Model.AddGe(E[k], A[i][j][k]);
                    }
                }
            }
            #endregion

            #region Constraint
            for (int i = 0; i < r_n + 1; i++)
            {
                for (int j = 1; j < r_n + 1; j++)
                {
                    for (int h = 0; h < r_n + 2; h++)
                    {
                        for (int k = 0; k < numtruck; k++)
                        {
                            if (i != j && j != h && h != i)
                                Model.AddGe(E[k], beta_ijg[i][j][h][k]);
                        }
                    }
                }
            }
            #endregion

            #region Constraint 
            for (int k = 0; k < numtruck; k++)
            {
                Model.AddEq(a_T[0][k], 0);
                Model.AddEq(a_D[0][k], 0);

            }
            #endregion
            #region Constraint 
            for (int k = 0; k < numtruck; k++)
            {
                for (int i = 0; i < r_n + 2; i++)
                {
                    Model.AddEq(X[i][k], Model.Max(a_D[i][k], a_T[i][k]));
                }
            }
            #endregion
            #region Constraint
            INumExpr[] expr222 = new INumExpr[1];
            for (int j = 1; j < r_n + 1; j++)
            {
                for (int k = 0; k < numtruck; k++)
                {
                    expr222[0] = beta_ijg[0][0][0][0];
                    expr222[0] = Model.Sum(expr222[0], Model.Prod(-1, beta_ijg[0][0][0][0]));
                    for (int i = 0; i < r_n + 1; i++)
                    {
                        expr222[0] = Model.Sum(expr222[0], A[i][j][k]);
                    }
                    Model.AddEq(p_ijk[0][j][k], expr222[0]);
                }
            }
            #endregion

            #region  objective function
            INumExpr[] ob = new INumExpr[numtruck];

            INumExpr[] obj2 = new INumExpr[1];
            INumExpr[] obj4 = new INumExpr[1];
            INumExpr[] obj6 = new INumExpr[1];
            INumExpr[] obj8 = new INumExpr[1];

            for (int k = 0; k < numtruck; k++)
            {
                obj2[0] = beta_ijg[0][1][2][0];
                obj2[0] = Model.Sum(obj2[0], Model.Prod(beta_ijg[0][1][2][0], -1));
                for (int j = 1; j < r_n + 1; j++)
                {
                    for (int h = 1; h < r_n + 2; h++)
                    {
                        for (int i = 0; i < r_n + 1; i++)
                        {
                            obj2[0] = Model.Sum(obj2[0], Model.Prod(beta_ijg[i][j][h][k], timeD[i, j] + timeD[j, h]));
                        }
                    }
                }

                obj4[0] = A[0][1][0];
                obj4[0] = Model.Sum(obj4[0], Model.Prod(A[0][1][0], -1));
                for (int j = 1; j < r_n + 2; j++)
                {
                    for (int i = 0; i < r_n + 1; i++)
                    {
                        obj4[0] = Model.Sum(obj4[0], Model.Prod(A[i][j][k], timeT[i, j, k]));
                    }
                }
                obj6[0] = A[0][1][0];
                obj6[0] = Model.Sum(obj6[0], Model.Prod(A[0][1][0], -1));
                obj6[0] = Model.Sum(obj6[0], X[r_n + 1][k]);
                obj8[0] = E[0];
                obj8[0] = Model.Sum(obj8[0], Model.Prod(E[0], -1));
                obj8[0] = Model.Sum(obj8[0], E[k]);
                ob[k] = Model.Sum(Model.Sum(Model.Sum(Model.Prod(ck - cw, obj4[0]), Model.Prod(cd, obj2[0])),
                    Model.Prod(cw, obj6[0])), Model.Prod(cg, obj8[0]));
            }

            INumExpr[] obj = new INumExpr[1];
            obj[0] = beta_ijg[0][1][2][0];
            obj[0] = Model.Sum(obj[0], Model.Prod(beta_ijg[0][1][2][0], -1));
            for (int k = 0; k < numtruck; k++)
            {
                obj[0] = Model.Sum(obj[0], ob[k]);
            }
            Model.AddMinimize(obj[0]);
            #endregion

            Model.ExportModel("TruckUAV.LP");

            Model.SetParam(Cplex.DoubleParam.TimeLimit, 10800);
            if (Model.Solve())
            {
                Console.WriteLine("cost=" + Math.Round(Model.ObjValue, 0));
                OBBBBJ = Math.Round(Model.ObjValue, 2);

                for (int k = 0; k < numtruck; k++)
                {
                    for (int i = 0; i < r_n + 2; i++)
                    {
                        for (int j = 0; j < r_n + 2; j++)
                        {
                            try
                            {
                                if (Model.GetValue(A[i][j][k]) > 0.5)
                                {
                                    Console.WriteLine("A[" + i.ToString() + "][" + j.ToString() + "][" + k.ToString() + "] = 1");
                                }
                            }
                            catch
                            {
                                continue;
                            }
                        }
                    }

                    Console.WriteLine();
                }
                Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                for (int i = 0; i < r_n + 1; i++)
                {
                    for (int j = 0; j < r_n + 1; j++)
                    {
                        for (int h = 0; h < r_n + 2; h++)
                        {
                            for (int d = 0; d < numtruck; d++)
                            {
                                try
                                {
                                    if (Model.GetValue(beta_ijg[i][j][h][d]) > 0.5)
                                    {
                                        Console.WriteLine("beta[" + i.ToString() + "][" + j.ToString() + "][" + h.ToString() + "][" + d.ToString() + "] = 1");
                                    }
                                }
                                catch
                                {
                                    continue;
                                }
                            }
                        }
                    }
                }
                Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                for (int k = 0; k < numtruck; k++)
                {
                    for (int i = 0; i < r_n + 2; i++)
                    {
                        for (int j = 0; j < r_n + 2; j++)
                        {

                            try
                            {
                                if (Model.GetValue(A[i][j][k]) >= 0.7)
                                {
                                    Console.WriteLine("a_T[" + i.ToString() + "][" + k.ToString() + "] = " + Math.Round(Model.GetValue(a_T[i][k]), 2));
                                    Console.WriteLine("a_T[" + j.ToString() + "][" + k.ToString() + "] = " + Math.Round(Model.GetValue(a_T[j][k]), 2));
                                    Console.WriteLine("timeT[" + i.ToString() + "][" + j.ToString() + "][" + k.ToString() + "] = " + timeT[i, j, k]);
                                    Console.WriteLine("\n");
                                }
                            }
                            catch
                            {
                                continue;
                            }
                        }
                    }
                }
                Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");

                for (int k = 0; k < numtruck; k++)

                {
                    for (int i = 0; i < r_n + 2; i++)

                    {
                        try
                        {
                            if (Model.GetValue(U[i][k]) != 0)
                            {
                                Console.WriteLine("U[" + i.ToString() + "][" + k.ToString() + "] = " + Model.GetValues(U[i])[k].ToString());
                            }
                        }
                        catch
                        {
                            continue;
                        }
                    }
                }
                Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                for (int k = 0; k < numtruck; k++)

                {
                    for (int i = 0; i < r_n + 2; i++)

                    {
                        for (int j = 1; j < r_n + 2; j++)

                        {
                            try
                            {
                                if (Model.GetValue(p_ijk[i][j][k]) != 0)
                                {
                                    Console.WriteLine("p_ijk[" + i.ToString() + "][" + j.ToString() + "][" + k.ToString() + "] = " + Model.GetValues(p_ijk[i][j])[k].ToString());
                                }
                            }
                            catch
                            {
                                continue;
                            }
                        }
                    }
                }
                Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                for (int d = 0; d < numtruck; d++)
                {
                    for (int i = 0; i < r_n + 1; i++)
                    {
                        for (int j = 1; j < r_n + 1; j++)
                        {
                            for (int h = 1; h < r_n + 2; h++)
                            {

                                try
                                {
                                    if (Model.GetValue(beta_ijg[i][j][h][d]) > 0.9)
                                    {
                                        Console.WriteLine("a_D[" + i.ToString() + "][" + d.ToString() + "] = " + Math.Round(Model.GetValue(a_D[i][d]), 2));
                                        Console.WriteLine("a_D[" + j.ToString() + "][" + d.ToString() + "] = " + Math.Round(Model.GetValue(a_D[j][d]), 2));
                                        Console.WriteLine("a_D[" + h.ToString() + "][" + d.ToString() + "] = " + Math.Round(Model.GetValue(a_D[h][d]), 2));
                                        Console.WriteLine("time_D[" + i.ToString() + "][" + j.ToString() + "][" + d.ToString() + "] = " + timeD[i, j]);
                                        Console.WriteLine("time_D[" + j.ToString() + "][" + h.ToString() + "][" + d.ToString() + "] = " + timeD[j, h]);
                                        Console.WriteLine("\n");
                                    }
                                }
                                catch
                                {
                                    continue;
                                }
                            }
                        }
                    }
                }
                Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                for (int k = 0; k < numtruck; k++)
                {
                    for (int i = 1; i < r_n + 2; i++)
                    {
                        try
                        {
                            if (Model.GetValue(X[i][k]) > 0)
                                Console.WriteLine("X[" + i + "][" + k + "]=" + Model.GetValue(X[i][k]));
                        }
                        catch
                        {
                            continue;
                        }
                    }
                }

                Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                for (int k = 0; k < numtruck; k++)
                {
                    for (int i = 0; i < r_n + 2; i++)
                    {
                        if (Model.GetValue(a_D[i][k]) != 0)

                            Console.WriteLine("aD[" + i + "][" + k + "]=" + Model.GetValue(a_D[i][k]));
                    }
                }
                Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                for (int k = 0; k < numtruck; k++)
                {
                    for (int i = 0; i < r_n + 2; i++)
                    {
                        if (Model.GetValue(a_T[i][k]) != 0)

                            Console.WriteLine("aT[" + i + "][" + k + "]=" + Model.GetValue(a_T[i][k]));
                    }
                }
                Model.End();
            }
            else
            {
                Model.End();
                Console.WriteLine();
                Console.WriteLine("cannot be solved");
            }
        }

        static void Main(string[] args)
        {
            Program prg = new Program();
            prg.input();
            double obj = 0;
            DateTime A = DateTime.Now;
            prg.solve(ref obj);
            Console.WriteLine("*************");
            DateTime final_time = DateTime.Now;
            TimeSpan CC = final_time - A;
            Console.WriteLine("Time={0}", CC);
            Console.WriteLine("obj={0}", obj);
            Console.ReadLine();
        }
    }
}
