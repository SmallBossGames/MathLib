package smallBossMathLib.explicitDifferentialEquations

import smallBossMathLib.differentialEquations.IFirstOrderIntegrator
import kotlin.math.*

class StabilityControlSecondOrderIntegrator(override val evaluations: Int, override val accuracy: Double)
    : RungeKuttaIntegratorBase(), IFirstOrderIntegrator
{
    override fun integrate(
        t0: Double,
        y0: DoubleArray,
        t: Double,
        outY: DoubleArray,
        equations: (t: Double, inY: DoubleArray, outY: DoubleArray) -> Unit
    ) {
        executeStepHandlers(t0, y0)

        y0.copyInto(outY)

        var fCurrentBuffer = DoubleArray(y0.size)
        var fLastBuffer = DoubleArray(y0.size)
        var currentEvaluationsCount = 0

        val yNextBuffer = DoubleArray(y0.size)

        val k1Buffer = DoubleArray(y0.size)
        val k2Buffer = DoubleArray(y0.size)
        val k3Buffer = DoubleArray(y0.size)

        var time = t0
        var step = t / evaluations

        equations(time, outY, fLastBuffer)

        while (time < t0 + t && isNextEvaluationAllow(currentEvaluationsCount)) {
            for (i in fLastBuffer.indices) {
                k1Buffer[i] = step * fLastBuffer[i]
                yNextBuffer[i] = outY[i] + 2.0 / 3.0 * k1Buffer[i]
            }

            equations(time, yNextBuffer, fCurrentBuffer)

            for (i in fCurrentBuffer.indices) {
                k2Buffer[i] = step * fCurrentBuffer[i]
                yNextBuffer[i] = outY[i] + 1.0 / 3.0 * k1Buffer[i] + 1.0 / 3.0 * k2Buffer[i]
            }

            equations(time, yNextBuffer, fCurrentBuffer)

            for (i in fCurrentBuffer.indices) {
                k3Buffer[i] = step * fCurrentBuffer[i]
                yNextBuffer[i] = outY[i] + 1.0 / 4.0 * k1Buffer[i] + 15.0 / 32.0 * k2Buffer[i] + 9.0 / 32.0 * k3Buffer[i]
            }
            
            val timeNext = step + time
            equations(timeNext, yNextBuffer, fCurrentBuffer)

            val q1 = findQ1(k1Buffer, k2Buffer, accuracy)
            if (q1 < 1.0) {
                step = min(t + t0 - time , q1 * step / 1.1)
                continue
            }

            val q2 = findQ2(fLastBuffer, fCurrentBuffer, step, accuracy)
            val r = findR(k1Buffer, k2Buffer, k3Buffer)


            if (q2 < 1.0 && r < 1) {
                step = min(t + t0 - time , q2 * step / 1.1)
                continue
            }

            step = if(q2 < 1){
                min(t + t0 - time , min(q1, q2)*step)
            } else{
                max(step, min(q1, min(q2, r))*step)
            }

            val temp = fLastBuffer
            fLastBuffer = fCurrentBuffer
            fCurrentBuffer = temp

            time = timeNext

            for (i in yNextBuffer.indices) {
                outY[i] = yNextBuffer[i]
            }

            executeStepHandlers(time, outY)

            currentEvaluationsCount++
        }
    }

    private fun findQ1(k1: DoubleArray, k2: DoubleArray, accuracy: Double) : Double{
        require(k1.size == k2.size)

        var norm = 0.0
        for (i in k1.indices){
            val d = (k2[i] - k1[i])
            norm += d*d
        }
        norm = sqrt(norm)

        //стр. 95
        return (((6.0 * 2.0/3.0 * accuracy) / (1.0 - 6.0 * 1.0/16.0)) / norm)
    }

    private fun findQ2(f1: DoubleArray, f2: DoubleArray, step:Double, accuracy: Double) : Double {
        require(f1.size == f2.size)

        var norm = 0.0
        for (i in f1.indices){
            val d = f1[i] - f2[i]
            norm += d*d
        }
        norm *= step

        return ((6.0*accuracy / (1.0 - 6.0 * 1.0 / 16.0)) / norm)
    }

    private fun findR(k1: DoubleArray, k2: DoubleArray, k3: DoubleArray) : Double {
        require(k1.size == k2.size && k2.size == k3.size)

        var max = Double.NEGATIVE_INFINITY
        for (i in k1.indices){
            val result = abs((k3[i] - k2[i]) / (k2[i] - k1[i]))
            if(result > max)
                max = result
        }

        return 2/max
    }
}