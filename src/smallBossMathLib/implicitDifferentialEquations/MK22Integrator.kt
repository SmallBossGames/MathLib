package smallBossMathLib.implicitDifferentialEquations

import smallBossMathLib.shared.*
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.min
import kotlin.math.sqrt

const val a = 0.29289321881
const val p1 = 5.0/4.0
const val p2 = 3.0/4.0
const val beta = 2.0/3.0
const val alpha = -4.0/3.0
//const val accuracyCoefficient = (6.0*a - 5.0)/(4.0 - 8.0*a)

const val rMin = 1e-14

/*const val a = 0.29289321881
const val p1 = 0.29289321881
const val p2 = 0.5/a
const val beta = 0.29289321881
const val alpha = -2*a
const val accuracyCoefficient = (6.0*a - 5.0)/(4.0 - 8.0*a)*/

const val v = 1

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

        if(y0.size != outY.size)
            throw IllegalArgumentException()

        executeStepHandlers(t0, y0)

        y0.copyInto(outY)

        val endTime = t0 + t

        val vectorBuffer1 = DoubleArray(y0.size)
        val vectorBuffer2 = DoubleArray(y0.size)
        val vectorBuffer3 = DoubleArray(y0.size)
        val vectorBuffer4 = DoubleArray(y0.size)

        val matrixBuffer = Matrix2D(y0.size)
        val jacobiMatrix = Matrix2D(y0.size)
        val dMatrix = Matrix2D(y0.size)

        var step = t / evaluations
        var time = t0
        var freezeSteps = 0
        var isNeedFindJacobi = true
        var currentEvaluationsCount = 0

        while (time < endTime && isNextEvaluationAllow(currentEvaluationsCount)){
            if(freezeSteps == 0){
                if(isNeedFindJacobi){
                    //findJacobiMatrix(outY, time, equations, jacobiMatrix)

                    outY.copyInto(vectorBuffer1)
                    for (i in vectorBuffer1.indices){
                        val r = max(rMin, sqrt(rMin)*abs(outY[i]))
                        vectorBuffer1[i] += r
                        equations(time, outY, vectorBuffer2)
                        equations(time, vectorBuffer1, vectorBuffer3)
                        for (j in vectorBuffer1.indices){
                            jacobiMatrix[j, i] = (vectorBuffer3[j] - vectorBuffer2[j]) / r
                        }
                        vectorBuffer1[i] = outY[i]
                    }
                }

                multipleMatrix2D(step * a, jacobiMatrix, dMatrix)
                for (i in dMatrix.indices){
                    dMatrix[i,i] = 1.0 - dMatrix[i,i]
                }
                dMatrix.makeLU()
            } else{
                isNeedFindJacobi = true
            }



            for (i in vectorBuffer1.indices){
                vectorBuffer1[i] = outY[i]
            }
            equations(time, vectorBuffer1, vectorBuffer2)
            for (i in vectorBuffer2.indices){
                vectorBuffer2[i] = step * vectorBuffer2[i]
            }
            dMatrix.solveLU(vectorBuffer2, vectorBuffer3)

            for (i in vectorBuffer1.indices){
                vectorBuffer1[i] = outY[i] + beta * vectorBuffer3[i]
            }
            equations(time, vectorBuffer1, vectorBuffer2)
            for (i in vectorBuffer2.indices){
                vectorBuffer2[i] = step * vectorBuffer2[i] + alpha * vectorBuffer3[i]
            }
            dMatrix.solveLU(vectorBuffer2, vectorBuffer4)

            /*var e1Norm = 0.0
            for (i in vectorBuffer1.indices){
                val temp = accuracyCoefficient * (vectorBuffer4[i] + 1.0/3.0*vectorBuffer3[i])
                vectorBuffer1[i] = temp
                val norm = abs(temp) / (abs(outY[i]) + v)
                if (norm > e1Norm)
                    e1Norm = norm
            }*/

            var e1Norm = 0.0
            for (i in vectorBuffer1.indices){
                val temp = (vectorBuffer4[i] + 1.0/3.0*vectorBuffer3[i])
                vectorBuffer1[i] = temp
                val norm = abs(temp) / (abs(outY[i]) + v)
                if (norm > e1Norm)
                    e1Norm = norm
            }

            val q1 = sqrt(8.0*accuracy / e1Norm)

            var e2Norm = 0.0
            val q2: Double
            if(q1 < 1.0) {
                dMatrix.inverseLU(matrixBuffer)
                matrixBuffer.multiply(vectorBuffer1, vectorBuffer2)
                for (i in vectorBuffer1.indices){
                    val norm = abs(vectorBuffer2[i]) / (abs(outY[i]) + v)
                    if (norm > e2Norm){
                        e2Norm = norm
                    }
                }
                q2 = sqrt(8.0*accuracy / e2Norm)
            } else {
                e2Norm = e1Norm
                q2 = q1
            }

            if(q2 < 1.0){
                val stepNew = step * q2
                step = normalizeStep(stepNew, time, endTime)

                if(freezeSteps==0){
                    isNeedFindJacobi = false
                } else {
                    freezeSteps = 0
                }
                continue
            }

            for (i in outY.indices){
                outY[i] = outY[i] + p1*vectorBuffer3[i] + p2*vectorBuffer4[i]
            }

            time += step
            freezeSteps++
            currentEvaluationsCount++

            executeStepHandlers(time, outY)

            val stepNew = min(q1, q2) * step

            if(!(freezeSteps < freezeJacobiSteps && stepNew < stepSizeCoefficient*step && e1Norm <= e2Norm)){
                freezeSteps = 0
                step = normalizeStep(stepNew, time, endTime)
            }
        }
    }

    private fun normalizeStep(step: Double, t: Double, endT: Double) =
        when {
            step < minStep -> {
                minStep
            }
            step > maxStep -> {
                maxStep
            }
            t + step > endT -> {
                endT - t
            }
            else -> {
                step
            }
        }
}