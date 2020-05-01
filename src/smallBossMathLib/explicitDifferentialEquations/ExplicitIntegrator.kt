package smallBossMathLib.explicitDifferentialEquations

import smallBossMathLib.shared.Integrator

abstract class ExplicitIntegrator : Integrator() {
    private val stepHandlers = HashSet<(
        time:Double,
        y:DoubleArray,
        state: IExplicitMethodStepState,
        statistic: IExplicitMethodStatistic
    ) -> Unit>()

    protected fun executeStepHandlers(
        time:Double,
        y:DoubleArray,
        state: IExplicitMethodStepState,
        statistic: IExplicitMethodStatistic
    ) {
        for (handler in stepHandlers)
            handler(time, y, state, statistic)
    }

    fun addStepHandler(handler: (
        time:Double,
        y:DoubleArray,
        state: IExplicitMethodStepState,
        statistic: IExplicitMethodStatistic
    ) -> Unit) {
        stepHandlers += handler
    }

    fun removeStepHandler(handler: (
        time:Double,
        y:DoubleArray,
        state: IExplicitMethodStepState,
        statistic: IExplicitMethodStatistic
    ) -> Unit) {
        stepHandlers -= handler
    }
}