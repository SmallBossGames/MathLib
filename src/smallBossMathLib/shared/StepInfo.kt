package smallBossMathLib.shared

class StepInfo(
    var time:Double,
    var yValue:DoubleArray,
    var isLowLimitReached:Boolean,
    var isHighLimitReached: Boolean){
    fun set(time:Double, yValue:DoubleArray, isLowLimitReached:Boolean, isHighLimitReached: Boolean){
        this.time = time
        this.yValue = yValue
        this.isLowLimitReached = isLowLimitReached
        this.isHighLimitReached = isHighLimitReached
    }
}