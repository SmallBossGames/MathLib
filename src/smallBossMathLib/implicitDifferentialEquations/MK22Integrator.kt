package smallBossMathLib.implicitDifferentialEquations

import smallBossMathLib.implicitDifferentialEquations.exceptions.ExceedingLimitEvaluationsException
import smallBossMathLib.implicitDifferentialEquations.exceptions.ExceedingLimitStepsException
import smallBossMathLib.shared.*
import kotlin.math.*

const val a = 0.29289321881
const val p1 = 0.29289321881
const val p2 = 0.5/a
const val beta = 0.29289321881
const val alpha = -2*a

class MK22Integrator (val startEvaluationCount: Int,
                      val freezeJacobiSteps: Int,
                      val stepSizeCoefficient: Double,
                      val accuracy: Double) : ImplicitIntegrator() {

    @Throws(ExceedingLimitStepsException::class, ExceedingLimitEvaluationsException::class)
    fun integrate(
        t0: Double,
        y0: DoubleArray,
        t: Double,
        rVector: DoubleArray,
        outY: DoubleArray,
        equations: (inY: DoubleArray, outF: DoubleArray) -> Unit) : IImplicitMethodStepInfo {

        if (y0.size != outY.size)
            throw IllegalArgumentException()

        y0.copyInto(outY)

        val endTime = t0 + t

        val vectorBuffer1 = DoubleArray(y0.size)
        val vectorBuffer2 = DoubleArray(y0.size)
        val k1 = DoubleArray(y0.size)
        val k2 = DoubleArray(y0.size)

        val matrixBuffer = Matrix2D(y0.size)
        val dMatrix = Matrix2D(y0.size)
        val jacobiMatrix = Matrix2D(y0.size)

        var step = t / startEvaluationCount
        var time = t0
        var freezeStepsCount = 0
        var isNeedFindJacobi = true

        val stepInfo = ImplicitMethodStepInfo(
            isLowStepSizeReached = isLowStepSizeReached(step),
            isHighStepSizeReached = isHighStepSizeReached(step),
            stepsCount = 0,
            evaluationsCount = 0,
            jacobiEvaluationsCount = 0,
            returnsCount = 0
        )

        executeStepHandlers(time, outY, stepInfo)

        while (time < endTime){
            step = normalizeStep(step, time, endTime)

            checkStepCount(stepInfo.stepsCount)

            if(freezeStepsCount == 0 && isNeedFindJacobi) {
                checkEvaluationCount(stepInfo.evaluationsCount)
                equations(outY, vectorBuffer2)
                stepInfo.evaluationsCount++
                outY.copyInto(vectorBuffer1)
                for (i in vectorBuffer1.indices){
                    val r = max(1e-14, 1e-7*abs(outY[i]))
                    val jacobiColumn = jacobiMatrix.columns[i]
                    vectorBuffer1[i] += r
                    checkEvaluationCount(stepInfo.evaluationsCount)
                    equations(vectorBuffer1, jacobiColumn)
                    stepInfo.evaluationsCount++
                    for (j in jacobiColumn.indices){
                        jacobiColumn[j] = (jacobiColumn[j] - vectorBuffer2[j]) / r
                    }
                    vectorBuffer1[i] = outY[i]
                }
                stepInfo.jacobiEvaluationsCount++
            } else {
                isNeedFindJacobi = true
            }

            for (i in dMatrix.indices)
                for(j in dMatrix.indices)
                    dMatrix[i, j] = (if (i == j) 1.0 else 0.0) - a*step*(jacobiMatrix[i,j])

            dMatrix.makeLU()

            for (i in vectorBuffer1.indices){
                vectorBuffer1[i] = outY[i]
            }
            checkEvaluationCount(stepInfo.evaluationsCount)
            equations(vectorBuffer1, vectorBuffer2)
            stepInfo.evaluationsCount++
            for (i in vectorBuffer2.indices){
                vectorBuffer2[i] = step*vectorBuffer2[i]
            }
            dMatrix.solveLU(vectorBuffer2, k1)
            for (i in vectorBuffer1.indices){
                vectorBuffer1[i] = outY[i] + beta*k1[i]
            }
            checkEvaluationCount(stepInfo.evaluationsCount)
            equations(vectorBuffer1, vectorBuffer2)
            stepInfo.evaluationsCount++
            for (i in vectorBuffer2.indices){
                vectorBuffer2[i] = step*vectorBuffer2[i] + alpha*k1[i]
            }
            dMatrix.solveLU(vectorBuffer2, k2)

            for (i in vectorBuffer1.indices){
                vectorBuffer1[i] = k2[i] + (2*a-1.0)*k1[i]
            }

            val e1 = zeroSafetyNorm(vectorBuffer1, outY, rVector)
            val q1 = sqrt(accuracy*(abs(a-2*a*a)/abs(a - 1.0/3.0)) / e1)

            val e2: Double
            val q2: Double
            if(q1 < 1.0){
                dMatrix.inverseLU(matrixBuffer)
                matrixBuffer.multiply(vectorBuffer1, vectorBuffer2)
                e2 = zeroSafetyNorm(vectorBuffer2, outY, rVector)
                q2 = sqrt(accuracy*(abs(a-2*a*a)/abs(a - 1.0/3.0)) / e2)
            } else {
                e2 = e1
                q2 = q1
            }

            if (q2 < 1.0 && !stepInfo.isLowStepSizeReached && !stepInfo.isHighStepSizeReached){
                if(freezeStepsCount == 0) {
                    isNeedFindJacobi = false
                } else {
                    freezeStepsCount = 0
                }

                step *= q2

                stepInfo.isLowStepSizeReached = isLowStepSizeReached(step)
                stepInfo.isHighStepSizeReached = isHighStepSizeReached(step)

                stepInfo.returnsCount++

                continue
            }

            for (i in outY.indices) {
                outY[i] = outY[i] + p1*k1[i] + p2*k2[i]
            }

            time += step
            stepInfo.stepsCount++

            executeStepHandlers(time, outY, stepInfo)

            val stepNew = min(q1, q2)*step

            freezeStepsCount++

            if(!(freezeStepsCount < freezeJacobiSteps && stepNew < step*stepSizeCoefficient && e1 <= e2)){
                stepInfo.isLowStepSizeReached = isLowStepSizeReached(step)
                stepInfo.isHighStepSizeReached = isHighStepSizeReached(step)

                step = stepNew
                freezeStepsCount = 0
            }
        }
        return stepInfo
    }
}