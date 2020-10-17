# using Plots
using Statistics
using DelimitedFiles
using Distributions

println("Welcome to this mini-script, testing your skills. You will be asked to hit 'Return' in 0.2seconds intervals and you will get a statistic of your precision and accuracy")
println("Press return to continue")
readline()
println("Choose the number of intervals you want to record (the more, the better). Type in the number and hit return afterwards")

try 
    global nmax = parse(Int64, readline())
catch e
    println("Please choose an Integer, start again!")
    exit()
end
println("You chose $nmax, press return to start; be ready to get going")
readline()

global n = 1;
times = []
# nmax = 100
while n < nmax
    println("Hit Return every 0.2s, $n / $nmax")
    time = @elapsed readline()
    println(time)
    push!(times, time)
    global n += 1
end

open("timerTest.dat", "w") do io
    writedlm(io, times)
end

println("You made it, 2s pause")
sleep(2)

println("Data saved as timerTest.dat")
println("Mean time ", mean(times), "s")
println("Standard deviation ", std(times), "s")
println("Hit return to continue")
d = Normal(mean(times), std(times))
x = range(minimum(times), maximum(times); length=500)

readline()

println("Do you want to see plots? Then type Y and hit return (might take time as Plots will be compiled)")
a = readline()
if a == "Y"
    using Plots
    p1 = plot(times, label="time series")
    p2 = histogram(times[5:end], normed=true, label="normed Histogram (omitting first 5)", bins=30)
    plot!(x, pdf.(d, x), label="Gaussian pdf")
    plot(p1, p2,  layout=(2, 1), legend=true)
    filename = "TimerTest.pdf"
    savefig(filename)
    println("Plot saved as $filename")
end
println("Hit return to finish")
readline()

