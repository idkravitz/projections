#include <octave/oct.h>
#include <octave/parse.h>
#include <octave/oct-convn.h>

#define EPSILON 1e-8

bool testIsNegative(double val)
{
        return val < -EPSILON;
}


DEFUN_DLD (SimProPartNative, args, nargout,
                "Native version of SimProPart")
{
        int nargin = args.length();
        Matrix X = args(0).matrix_value();
        Matrix X1, X2;
        octave_function *fcn = args(1).function_value();

        std::string octave_rand = "rand";
        std::string octave_norm = "norm";
        std::string simPro = "SimPro";

        do {
                octave_value_list randArgs;
                randArgs(0) = 1;
                randArgs(1) = X.columns();
                octave_value_list _w = feval(octave_rand, randArgs, 1);
                NDArray w = _w(0).array_value();

                octave_value_list normArgs;
                normArgs(0) = w;
                octave_value_list _normW = feval(octave_norm, normArgs, 1);
                double normW = _normW(0).double_value();

                w = w / normW;

                octave_value_list splitFuncArgs;
                splitFuncArgs(0) = X;
                splitFuncArgs(1) = X * ColumnVector(w);
                octave_value_list split = feval(fcn, splitFuncArgs, 2);

                X1 = split(0).matrix_value();
                X2 = split(1).matrix_value();

        } while(X1.columns() == 0 || X2.columns() == 0);

        octave_value_list simProArgs;
        simProArgs(1) = 1e5;
        simProArgs(2) = EPSILON;
        RowVector verboseLevel(2, -1);
        simProArgs(3) = verboseLevel;
        simProArgs(4) = NDArray();
        simProArgs(5) = NDArray();

        simProArgs(0) = X1;
        octave_value_list _x1 = feval(simPro, simProArgs, 1);
        ColumnVector x1 = _x1(0).column_vector_value();

        simProArgs(0) = X2;
        octave_value_list _x2 = feval(simPro, simProArgs, 1);
        ColumnVector x2 = _x2(0).column_vector_value();

        simProArgs(0) = Matrix(x1).append(x2);
        octave_value_list _z = feval(simPro, simProArgs, 1);
        RowVector z = _z(0).array_value();

        double zsmsq = NDArray(z).sumsq()(0);

        int it = 1;
        octave_value_list result;
        while ((z * X - zsmsq).test_any(testIsNegative) && it < 1000) {
                octave_value_list splitFuncArgs;
                ColumnVector zc = z.transpose();
                splitFuncArgs(0) = X;
                splitFuncArgs(1) = zc;
                octave_value_list split = feval(fcn, splitFuncArgs, 2);
                X1 = split(0).matrix_value();
                X2 = split(1).matrix_value();

                simProArgs(0) = X1.append(zc);
                octave_value_list _x1 = feval(simPro, simProArgs, 1);
                ColumnVector x1 = _x1(0).column_vector_value();

                simProArgs(0) = X2.append(zc);
                octave_value_list _x2 = feval(simPro, simProArgs, 1);
                ColumnVector x2 = _x2(0).column_vector_value();

                simProArgs(0) = Matrix(x1).append(x2);
                result = feval(simPro, simProArgs, nargout);
                z = result(0).array_value();
                zsmsq = NDArray(z).sumsq()(0);
                it++;
        }

        return result;
}
