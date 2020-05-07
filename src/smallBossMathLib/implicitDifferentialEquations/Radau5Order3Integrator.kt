package smallBossMathLib.implicitDifferentialEquations

import smallBossMathLib.implicitDifferentialEquations.exceptions.ExceedingLimitEvaluationsException
import smallBossMathLib.implicitDifferentialEquations.exceptions.ExceedingLimitStepsException
import smallBossMathLib.shared.NewtonRaphsonSolver
import smallBossMathLib.shared.StationaryODE

class Radau5Order3Integrator(val accuracy: Double) : ImplicitIntegrator() {
    @Throws(ExceedingLimitStepsException::class, ExceedingLimitEvaluationsException::class)
    fun integrate(
        startTime: Double,
        endTime: Double,
        defaultStepSize: Double,
        y0: DoubleArray,
        rVector: DoubleArray,
        outY: DoubleArray,
        equations: StationaryODE
    ) : IImplicitMethodStatistic {
        if (y0.size != outY.size)
            throw IllegalArgumentException()

        y0.copyInto(outY)

        var fCurrentBuffer = DoubleArray(y0.size)
        var fNextBuffer = DoubleArray(y0.size)

        val nextY = DoubleArray(y0.size)
        val vectorBuffer1 = DoubleArray(y0.size)

        var step = defaultStepSize
        var time = startTime

        val newtonSolver = NewtonRaphsonSolver(y0.size, accuracy, rVector)

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

        equations(y0, fCurrentBuffer)

        return statistic
    }
}