package smallBossMathLib.shared

import smallBossMathLib.implicitDifferentialEquations.exceptions.ExceedingLimitEvaluationsException
import smallBossMathLib.implicitDifferentialEquations.exceptions.ExceedingLimitStepsException

abstract class IntegratorBase {
    var isStepCountLimitEnabled: Boolean = false
        private set

    var isEvaluationsCountLimitEnabled: Boolean = false
        private set 

    var maxStepCount: Int = 0
        private set
    
    var maxEvaluationsCount: Int = 0
        private set

    private val stepHandlers = HashSet<(t: Double, y: DoubleArray) -> Unit>()

    fun addStepHandler(handler: (t: Double, y: DoubleArray) -> Unit) {
        stepHandlers += handler
    }

    fun removeStepHandler(handler: (t: Double, y: DoubleArray) -> Unit) {
        stepHandlers -= handler
    }

    fun clearStepHandlers() {
        stepHandlers.clear()
    }

    protected fun executeStepHandlers(t: Double, y: DoubleArray) {
        for (handler in stepHandlers)
            handler(t, y)
    }

    //Evaluation control
    fun enableStepCountLimit(maxSteps: Int){
        isStepCountLimitEnabled = true
        maxStepCount = maxSteps
    }

    fun disableStepCountLimit(){
        isStepCountLimitEnabled = false
        maxStepCount = 0
    }

    @Throws(ExceedingLimitStepsException::class)
    fun checkStepCount(steps: Int){
        if(isStepCountLimitEnabled && steps >= maxStepCount){
            throw ExceedingLimitStepsException()
        }
    }

    //Right part evaluations control

    fun enableEvaluationsCountLimit(maxEvaluations: Int){
        isStepCountLimitEnabled = true
        maxStepCount = maxEvaluations
    }

    fun disableEvaluationsCountLimit(){
        isStepCountLimitEnabled = false
        maxStepCount = 0
    }

    @Throws(ExceedingLimitEvaluationsException::class)
    fun checkEvaluationCount(evaluations: Int){
        if(isEvaluationsCountLimitEnabled && evaluations >= maxEvaluationsCount) {
            throw ExceedingLimitEvaluationsException()
        }
    }

    /*protected fun isNextStepAllowed(steps: Int, evaluations: Int) : Boolean =
        (!isStepCountLimitEnabled || steps < maxStepCount) &&
                (!isEvaluationsCountLimitEnabled || evaluations < maxEvaluationsCount)*/
}