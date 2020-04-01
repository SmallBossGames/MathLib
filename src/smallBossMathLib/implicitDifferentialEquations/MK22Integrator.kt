package smallBossMathLib.implicitDifferentialEquations

import smallBossMathLib.shared.*
import kotlin.math.*

/*const val a = 0.29289321881
const val p1 = 5.0/4.0
const val p2 = 3.0/4.0
const val beta = 2.0/3.0
const val alpha = -4.0/3.0*/

const val a = 0.29289321881
const val p1 = 0.29289321881
const val p2 = 0.5/a
const val beta = 0.29289321881
const val alpha = -2*a

const val v = 1e-7

class MK22Integrator (val evaluations: Int,
                      val freezeJacobiSteps: Int,
                      val stepSizeCoefficient: Double,
                      val accuracy: Double,
                      val minStep: Double,
                      val maxStep: Double) : IntegratorBase() {
    fun integrate(
        t0: Double,
        y0: DoubleArray,
        t: Double,
        outY: DoubleArray,
        equations: (t: Double, inY: DoubleArray, outF: DoubleArray) -> Unit) {

        if (y0.size != outY.size)
            throw IllegalArgumentException()

        executeStepHandlers(t0, y0)

        y0.copyInto(outY)

        val endTime = t0 + t

        val vectorBuffer1 = DoubleArray(y0.size)
        val vectorBuffer2 = DoubleArray(y0.size)
        val k1 = DoubleArray(y0.size)
        val k2 = DoubleArray(y0.size)

        val matrixBuffer = Matrix2D(y0.size)
        val dMatrix = Matrix2D(y0.size)
        val jacobiMatrix = Matrix2D(y0.size)
        val guessMatrix = Matrix2D(y0.size)

        var step = t / evaluations
        var time = t0
        var freezeStepsCount = 0
        var isNeedFindJacobi = true
        var currentEvaluationsCount = 0


        while (time < endTime && isNextEvaluationAllow(currentEvaluationsCount)){
            if(freezeStepsCount == 0 && isNeedFindJacobi) {
                    equations(time, outY, vectorBuffer2)
                    outY.copyInto(vectorBuffer1)
                    for (i in vectorBuffer1.indices){
                        val r = max(1e-14, 1e-7*abs(outY[i]))
                        val jacobiColumn = jacobiMatrix.columns[i]
                        val guessColumn = guessMatrix.columns[i]
                        vectorBuffer1[i] += r
                        equations(time, vectorBuffer1, jacobiColumn)
                        vectorBuffer1[i] += r
                        equations(time, vectorBuffer1, guessColumn)
                        for (j in vectorBuffer1.indices){
                            val jacobiValue1 = (jacobiColumn[j] - vectorBuffer2[j]) / r
                            val jacobiValue2 = (guessColumn[j] - jacobiColumn[j]) / r
                            val guessValue = 0.5 * (r / step) * (jacobiValue2 - jacobiValue1) / r
                            jacobiColumn[j] = jacobiValue1
                            guessColumn[j] = guessValue
                        }
                        vectorBuffer1[i] = outY[i]
                    }
            } else {
                isNeedFindJacobi = true
            }

            for (i in dMatrix.indices)
                for(j in dMatrix.indices)
                    dMatrix[i, j] = (if (i == j) 1.0 else 0.0) - a*step*(jacobiMatrix[i,j] + step * guessMatrix[i,j])

            dMatrix.makeLU()

            for (i in vectorBuffer1.indices){
                vectorBuffer1[i] = outY[i]
            }
            equations(time, vectorBuffer1, vectorBuffer2)
            for (i in vectorBuffer2.indices){
                vectorBuffer2[i] = step*vectorBuffer2[i]
            }
            dMatrix.solveLU(vectorBuffer2, k1)

            for (i in vectorBuffer1.indices){
                vectorBuffer1[i] = outY[i] + beta*k1[i]
            }
            equations(time, vectorBuffer1, vectorBuffer2)
            for (i in vectorBuffer2.indices){
                vectorBuffer2[i] = step*vectorBuffer2[i] + alpha*k1[i]
            }
            dMatrix.solveLU(vectorBuffer2, k2)

            var e1 = 0.0
            for (i in dMatrix.indices){
                val temp = k2[i] + (2*a-1.0)*k1[i]
                vectorBuffer1[i] = temp
                e1 = max(abs(temp) / (abs(outY[i]) + v), e1)
            }

            val q1 = sqrt(accuracy / e1)

            var e2 = 0.0
            val q2: Double
            if(q1 < 1.0){
                dMatrix.inverseLU(matrixBuffer)
                matrixBuffer.multiply(vectorBuffer1, vectorBuffer2)
                for (i in dMatrix.indices)
                    e2 = max(abs(vectorBuffer2[i]) / (abs(outY[i]) + v), e2)
                q2 = sqrt(accuracy / e2)
            } else {
                e2 = e1
                q2 = q1
            }

            if (q2 < 1.0){
                if(freezeStepsCount == 0) {
                    isNeedFindJacobi = false
                } else {
                    freezeStepsCount = 0
                }

                step = normalizeStep(step*q2, time, endTime)

                continue
            }

            for (i in outY.indices) {
                outY[i] = outY[i] + p1*k1[i] + p2*k2[i]
            }

            time += step
            currentEvaluationsCount++

            executeStepHandlers(time, outY)

            val hNew = min(q1, q2)*step

            freezeStepsCount++

            if(!(freezeStepsCount < freezeJacobiSteps && hNew < step*stepSizeCoefficient && e1 <= e2)){
                step = normalizeStep(hNew, time, endTime)
                freezeStepsCount = 0
            }
        }
    }

    private fun normalizeStep(step: Double, t: Double, endT: Double) : Double{
        val maxByLimit = if (t + step > endT) endT - t else step
        return when {
            maxByLimit < minStep -> minStep
            maxByLimit > maxStep -> maxStep
            else -> maxByLimit
        }
    }
}