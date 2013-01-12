#include <pthread.h>
#include <octave/oct.h>
#include <octave/parse.h>
#include <octave/oct-convn.h>

#define EPSILON 1e-8

bool testIsNegative(double val)
{
        return val < -EPSILON;
}

std::string simPro = "SimPro";
ColumnVector zc;
Matrix tasksMatrices[2];
ColumnVector retVals[2];
octave_value_list tasksArgs[2];

void *computationalTask(void *taskid) {
        int id = (int) taskid;
        octave_value_list &args = tasksArgs[id];
        args(0) = tasksMatrices[id].append(zc);
        octave_value_list _retVal = feval(simPro, args, 1);
        retVals[id] = _retVal(0).column_vector_value();

        pthread_exit(NULL);
}


DEFUN_DLD (SimProPartParallelNative, args, nargout,
                "Parallel Native version of SimProPart")
{
        int nargin = args.length();
        Matrix X = args(0).matrix_value();
//        Matrix &tasksMatrices[0] = tasksMatrices[0], &tasksMatrices[1] = tasksMatrices[1];
        octave_function *fcn = args(1).function_value();

        std::string octave_rand = "rand";
        std::string octave_norm = "norm";

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

                tasksMatrices[0] = split(0).matrix_value();
                tasksMatrices[1] = split(1).matrix_value();

        } while(tasksMatrices[0].columns() == 0 || tasksMatrices[1].columns() == 0);

        octave_value_list simProArgs;
        simProArgs(1) = 1e5;
        simProArgs(2) = EPSILON;
        RowVector verboseLevel(2, 1);
        simProArgs(3) = verboseLevel;
        simProArgs(4) = NDArray();
        simProArgs(5) = NDArray();

        simProArgs(0) = tasksMatrices[0];
        octave_value_list _x1 = feval(simPro, simProArgs, 1);
        ColumnVector x1 = _x1(0).column_vector_value();

        simProArgs(0) = tasksMatrices[1];
        octave_value_list _x2 = feval(simPro, simProArgs, 1);
        ColumnVector x2 = _x2(0).column_vector_value();

        simProArgs(0) = Matrix(x1).append(x2);
        octave_value_list _z = feval(simPro, simProArgs, 1);
        RowVector z = _z(0).array_value();

        double zsmsq = NDArray(z).sumsq()(0);
        for (int i = 1; i < simProArgs.length(); ++i) {
            tasksArgs[0](i) = simProArgs(i);
            tasksArgs[1](i) = simProArgs(i);
        }
        int it = 1;
        octave_value_list result;
        pthread_t thread_1, thread_2;
        while ((z * X - zsmsq).test_any(testIsNegative) && it < 1000) {
                octave_stdout << "loop start\n";
                octave_value_list splitFuncArgs;
                zc = z.transpose();
                splitFuncArgs(0) = X;
                splitFuncArgs(1) = zc;
                octave_value_list split = feval(fcn, splitFuncArgs, 2);
                tasksMatrices[0] = split(0).matrix_value();
                tasksMatrices[1] = split(1).matrix_value();

//                computationalTask((void *) 0);
//                computationalTask((void *) 1);

                pthread_create(&thread_1, NULL, computationalTask, (void *) 0);
                pthread_create(&thread_2, NULL, computationalTask, (void *) 1);

                pthread_join(thread_1, NULL);
                pthread_join(thread_2, NULL);
                octave_stdout << "we joined?\n";

                simProArgs(0) = Matrix(retVals[0]).append(retVals[1]);
                result = feval(simPro, simProArgs, nargout);
                z = result(0).array_value();
                zsmsq = NDArray(z).sumsq()(0);
                it++;
        }

        return result;
}
