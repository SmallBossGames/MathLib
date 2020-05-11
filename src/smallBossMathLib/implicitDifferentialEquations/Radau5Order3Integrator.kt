package smallBossMathLib.implicitDifferentialEquations

import smallBossMathLib.implicitDifferentialEquations.exceptions.ExceedingLimitEvaluationsException
import smallBossMathLib.implicitDifferentialEquations.exceptions.ExceedingLimitStepsException
import smallBossMathLib.shared.*
import kotlin.math.max
import kotlin.math.pow

const val uround = 1e-16

class Radau5Order3Integrator(
    val relativeTolerance: Double,
    val absoluteTolerance: Double,
    val maxNewtonIterations: Int
) : ImplicitIntegrator() {
    private companion object{
        val aMatrix : Matrix2D = Matrix2D(3)
        val bVector = doubleArrayOf(0.3764030627, 0.51248582618, 0.11111111111)
        val cVector = doubleArrayOf(0.15505102572, 0.64494897427, 1.0)
        val dVector = doubleArrayOf(0.0, 0.0, 1.0)
        val aMatrixEigenvalue = 0.274889

        init {
            aMatrix[0,0] = 0.19681547722
            aMatrix[0,1] = -0.06553542585
            aMatrix[0,2] = 0.02377097434
            aMatrix[1,0] = 0.39442431473
            aMatrix[1,1] = 0.29207341166
            aMatrix[1,2] = -0.04154875212
            aMatrix[2,0] = 0.3764030627
            aMatrix[2,1] = 0.51248582618
            aMatrix[2,2] = 0.11111111111
        }
    }

    @Throws(ExceedingLimitStepsException::class, ExceedingLimitEvaluationsException::class)
    fun integrate(
        startTime: Double,
        endTime: Double,
        defaultStepSize: Double,
        y0: DoubleArray,
        rVector: DoubleArray,
        outY: DoubleArray,
        equations: NonStationaryODE
    ) : IImplicitMethodStatistic {
        if (y0.size != outY.size)
            throw IllegalArgumentException()

        y0.copyInto(outY)

        val jacobiMatrix = Matrix2D(y0.size)
        val gMatrix = Matrix2D(aMatrix.size * jacobiMatrix.size)


        val identityMatrix = Matrix2D.createUnitMatrix2D(y0.size)
        val aIdentityMatrixBuffer = Matrix2D(gMatrix.size)

        aMatrix.kroneckerMultiply(identityMatrix, aIdentityMatrixBuffer)

        val fBuffer = DoubleArray(y0.size)
        val vectorBuffer = DoubleArray(y0.size)

        val zVector = DoubleArray(gMatrix.size)
        val vectorExtendedBuffer1 = DoubleArray(gMatrix.size)
        val vectorExtendedBuffer2 = DoubleArray(gMatrix.size)

        var step = defaultStepSize
        var time = startTime

        val newtonSolver = NewtonRaphsonSolver(y0.size, relativeTolerance, rVector)
        val jacobiSolver = JacobiMatrixSolver(y0.size)

        val statistic = ImplicitMethodStatistic(
            stepsCount = 0,
            evaluationsCount = 0,
            jacobiEvaluationsCount = 0,
            returnsCount = 0
        )

        val state = ImplicitMethodStepState(
            jacobiMatrix = newtonSolver.jacobiMatrix,
            isLowStepSizeReached = isLowStepSizeReached(step),
            isHighStepSizeReached = isHighStepSizeReached(step),
            freezeJacobiStepsCount = 0
        )

        var faccon = 1.0

        while (time < endTime){
            step = normalizeStep(step, time, endTime)

            jacobiSolver.solve(time, outY, jacobiMatrix, equations)
            aMatrix.kroneckerMultiply(jacobiMatrix, gMatrix)

            for (i in gMatrix.indices){
                for (j in gMatrix.indices){
                    gMatrix[i, j] = if (i == j) 1.0 else 0.0 - step * gMatrix[i, j]
                }
            }

            gMatrix.makeLU()

            for (i in zVector.indices){
                zVector[i] = 0.0
            }

            var newtonIterationsCount = 0
            var deltaZNormOld = 0.0

            while (true){
                for (i in aMatrix.indices){
                    for (j in vectorBuffer.indices){
                        vectorBuffer[j] = outY[j] + zVector[i * vectorBuffer.size + j]
                    }
                    equations(cVector[i] + time, vectorBuffer, fBuffer)
                    for (j in fBuffer.indices){
                        vectorExtendedBuffer1[i * vectorBuffer.size + j] = step * fBuffer[j]
                    }
                }

                aIdentityMatrixBuffer.multiply(vectorExtendedBuffer1, vectorExtendedBuffer2)

                for (i in vectorExtendedBuffer2.indices){
                    vectorExtendedBuffer2[i] = vectorExtendedBuffer2[i] - zVector[i]
                }

                gMatrix.solveLU(vectorExtendedBuffer2, vectorExtendedBuffer1)

                for (i in zVector.indices){
                    zVector[i] += vectorExtendedBuffer1[i]
                }

                val deltaZNorm = zeroSafetyNorm(vectorExtendedBuffer1, zVector, uround)

                faccon = if(newtonIterationsCount > 1){
                    val theta = deltaZNorm / deltaZNormOld
                    theta / (1.0 - theta)
                } else {
                    max(faccon, uround).pow(0.8)
                }

                newtonIterationsCount++

                if(faccon * deltaZNorm <= relativeTolerance * 5e-2)
                    break

                deltaZNormOld = deltaZNorm
            }

            for (i in zVector.indices){
                outY[i % outY.size] += dVector[i % dVector.size] * zVector[i]
            }

            time += step

            executeStepHandlers(time, outY, state, statistic)
        }
        return statistic
    }
}