
using PyPlot
using Statistics
using Printf
using DelimitedFiles

function initialize(Length, initial_state="random")

    if initial_state == "random"
        Config = rand(Bool, Length,Length)
    elseif initial_state == "ordered"
        Config = ones(Bool, Length,Length)
    end

    return Config
end

function stream(Config, Temperature, num_sample, num_sample_burn, period_sample)
    Length = size(Config)[1]

    vec_magnetization = zeros(num_sample)
    vec_energy = zeros(num_sample)

    for time = 1:num_sample_burn
        for iteration = 1:Length^2
            x = mod(rand(Int8),Length) + 1
            y = mod(rand(Int8),Length) + 1

            spin_site = Config[x,y]
            E = 4
            if Config[mod(x,Length)+1,y] == spin_site
                E -= 2
            end
            if Config[mod(x-2,Length)+1,y] == spin_site
                E -= 2
            end
            if Config[x,mod(y,Length)+1] == spin_site
                E -= 2
            end
            if Config[x,mod(y-2,Length)+1] == spin_site
                E -= 2
            end
            
            if E > 0
                Config[x,y] = !Config[x,y] # Flip
            elseif rand(Float64) < exp(2*E/Temperature)
                Config[x,y] = !Config[x,y] # Flip
            end
        end
    end
    magnetization = 0
    magnetization2 = 0
    for time = 1:num_sample
        for iteration = 1:Length^2*period_sample
            x = mod(rand(Int8),Length) + 1
            y = mod(rand(Int8),Length) + 1

            spin_site = Config[x,y]
            E = 4
            if Config[mod(x,Length)+1,y] == spin_site
                E -= 2
            end
            if Config[mod(x-2,Length)+1,y] == spin_site
                E -= 2
            end
            if Config[x,mod(y,Length)+1] == spin_site
                E -= 2
            end
            if Config[x,mod(y-2,Length)+1] == spin_site
                E -= 2
            end
            
            if E > 0
                Config[x,y] = !Config[x,y] # Flip
            elseif rand(Float64) < exp(2*E/Temperature)
                Config[x,y] = !Config[x,y] # Flip
            end
        end

        vec_magnetization[time] = -Length^2 + 2*sum(Config)
        vec_energy[time] = 4*Length^2 - 2*sum((Config.==circshift(Config, (-1,0))) + (Config.==circshift(Config, (1,0))) + (Config.==circshift(Config, (0,-1))) + (Config.==circshift(Config, (0,1))))
    end
    
    return vec_magnetization, vec_energy
end

function stream_tw(Config, Temperature, num_sample, num_sample_burn, period_sample)
    Length = size(Config)[1]

    vec_magnetization = zeros(num_sample)
    vec_energy = zeros(num_sample)

    for time = 1:num_sample_burn
        for x = 1:Length
            for y = 1:Length
                spin_site = Config[x,y]
                E = 4
                if Config[mod(x,Length)+1,y] == spin_site
                    E -= 2
                end
                if Config[mod(x-2,Length)+1,y] == spin_site
                    E -= 2
                end
                if Config[x,mod(y,Length)+1] == spin_site
                    E -= 2
                end
                if Config[x,mod(y-2,Length)+1] == spin_site
                    E -= 2
                end

                if E > 0
                    Config[x,y] = !Config[x,y] # Flip
                elseif rand(Float64) < exp(2*E/Temperature)
                    Config[x,y] = !Config[x,y] # Flip
                end
            end
        end
    end
    magnetization = 0
    magnetization2 = 0
    for time = 1:num_sample
        for iteration = 1:period_sample
            for x = 1:Length
                for y = 1:Length
                    spin_site = Config[x,y]
                    E = 4
                    if Config[mod(x,Length)+1,y] == spin_site
                        E -= 2
                    end
                    if Config[mod(x-2,Length)+1,y] == spin_site
                        E -= 2
                    end
                    if Config[x,mod(y,Length)+1] == spin_site
                        E -= 2
                    end
                    if Config[x,mod(y-2,Length)+1] == spin_site
                        E -= 2
                    end

                    if E > 0
                        Config[x,y] = !Config[x,y] # Flip
                    elseif rand(Float64) < exp(2*E/Temperature)
                        Config[x,y] = !Config[x,y] # Flip
                    end
                end
            end
        end

        vec_magnetization[time] = -Length^2 + 2*sum(Config)
        vec_energy[time] = 4*Length^2 - 2*sum((Config.==circshift(Config, (-1,0))) + (Config.==circshift(Config, (1,0))) + (Config.==circshift(Config, (0,-1))) + (Config.==circshift(Config, (0,1))))
    end
    
    return vec_magnetization, vec_energy
end

function stream_glauber(Config, Temperature, num_sample, num_sample_burn, period_sample)
    Length = size(Config)[1]

    vec_magnetization = zeros(num_sample)
    vec_energy = zeros(num_sample)

    for time = 1:num_sample_burn
        for iteration = 1:Length^2
            x = mod(rand(Int8),Length) + 1
            y = mod(rand(Int8),Length) + 1

            x = mod(rand(Int8),Length) + 1
            y = mod(rand(Int8),Length) + 1

            field =  -4 + 2*(Config[mod(x,Length)+1,y] + Config[mod(x-2,Length)+1,y] + Config[x,mod(y,Length)+1] + Config[x,mod(y-2,Length)+1])
            
            if rand(Float64) < 1 /(1 + exp(2*field/Temperature))
                Config[x,y] = 0
            else
                Config[x,y] = 1
            end
        end
    end
    magnetization = 0
    magnetization2 = 0
    for time = 1:num_sample
        for iteration = 1:Length^2*period_sample
            x = mod(rand(Int8),Length) + 1
            y = mod(rand(Int8),Length) + 1

            field =  -4 + 2*(Config[mod(x,Length)+1,y] + Config[mod(x-2,Length)+1,y] + Config[x,mod(y,Length)+1] + Config[x,mod(y-2,Length)+1])
            
            if rand(Float64) < 1 /(1 + exp(2*field/Temperature))
                Config[x,y] = 0
            else
                Config[x,y] = 1
            end
        end

        vec_magnetization[time] = -Length^2 + 2*sum(Config)
        vec_energy[time] = 4*Length^2 - 2*sum((Config.==circshift(Config, (-1,0))) + (Config.==circshift(Config, (1,0))) + (Config.==circshift(Config, (0,-1))) + (Config.==circshift(Config, (0,1))))
    end
    
    return vec_magnetization, vec_energy
end

function stream_glauber_tw(Config, Temperature, num_sample, num_sample_burn, period_sample)
    Length = size(Config)[1]

    vec_magnetization = zeros(num_sample)
    vec_energy = zeros(num_sample)

    for time = 1:num_sample_burn
        for x = 1:Length
            for y = 1:Length
                x = mod(rand(Int8),Length) + 1
                y = mod(rand(Int8),Length) + 1

                field =  -4 + 2*(Config[mod(x,Length)+1,y] + Config[mod(x-2,Length)+1,y] + Config[x,mod(y,Length)+1] + Config[x,mod(y-2,Length)+1])

                if rand(Float64) < 1 /(1 + exp(2*field/Temperature))
                    Config[x,y] = 0
                else
                    Config[x,y] = 1
                end
            end
        end
    end
    magnetization = 0
    magnetization2 = 0
    for time = 1:num_sample
        for iteration = 1:period_sample
            for x = 1:Length
                for y = 1:Length
                    x = mod(rand(Int8),Length) + 1
                    y = mod(rand(Int8),Length) + 1

                    field =  -4 + 2*(Config[mod(x,Length)+1,y] + Config[mod(x-2,Length)+1,y] + Config[x,mod(y,Length)+1] + Config[x,mod(y-2,Length)+1])

                    if rand(Float64) < 1 /(1 + exp(2*field/Temperature))
                        Config[x,y] = 0
                    else
                        Config[x,y] = 1
                    end
                end
            end
        end

        vec_magnetization[time] = -Length^2 + 2*sum(Config)
        vec_energy[time] = 4*Length^2 - 2*sum((Config.==circshift(Config, (-1,0))) + (Config.==circshift(Config, (1,0))) + (Config.==circshift(Config, (0,-1))) + (Config.==circshift(Config, (0,1))))
    end
    
    return vec_magnetization, vec_energy
end

function plot_config(Config, Temperature, num_sample, num_sample_burn, period_sample)
    Length = size(Config)[1]

    for time = 1:num_sample_burn
        for iteration = 1:Length^2
            x = mod(rand(Int8),Length) + 1
            y = mod(rand(Int8),Length) + 1

            spin_site = Config[x,y]
            E = 4
            if Config[mod(x,Length)+1,y] == spin_site
                E -= 2
            end
            if Config[mod(x-2,Length)+1,y] == spin_site
                E -= 2
            end
            if Config[x,mod(y,Length)+1] == spin_site
                E -= 2
            end
            if Config[x,mod(y-2,Length)+1] == spin_site
                E -= 2
            end
            
            if E > 0
                Config[x,y] = !Config[x,y] # Flip
            elseif rand(Float64) < exp(2*E/Temperature)
                Config[x,y] = !Config[x,y] # Flip
            end
        end
    end
    for time = 1:num_sample
        for iteration = 1:Length^2*period_sample
            x = mod(rand(Int8),Length) + 1
            y = mod(rand(Int8),Length) + 1

            spin_site = Config[x,y]
            E = 4
            if Config[mod(x,Length)+1,y] == spin_site
                E -= 2
            end
            if Config[mod(x-2,Length)+1,y] == spin_site
                E -= 2
            end
            if Config[x,mod(y,Length)+1] == spin_site
                E -= 2
            end
            if Config[x,mod(y-2,Length)+1] == spin_site
                E -= 2
            end
            
            if E > 0
                Config[x,y] = !Config[x,y] # Flip
            elseif rand(Float64) < exp(2*E/Temperature)
                Config[x,y] = !Config[x,y] # Flip
            end
        end
    end
    
    figure()
    imshow(Config, cmap="gray")
end

Length = 256
Config = initialize(Length, "random")

num_sample = 1
num_sample_burn = 4000
period_sample = 1

vec_Temperature = [1.5, 2.27, 3]

for index_Temperature = 1:length(vec_Temperature)
    Temperature = vec_Temperature[index_Temperature]
    print(Temperature, "  ")

    plot_config(Config, Temperature, num_sample, num_sample_burn, period_sample)
end    

Length = 16
Config = initialize(Length, "random")

num_sample = 1000
num_sample_burn = 100

vec_Temperature = vcat(1.5:0.1:2.2, 2.21:0.01:2.49, 2.5:0.1:3) 
data_magnetization = zeros(length(vec_Temperature))
data_chi = zeros(length(vec_Temperature))
data_energy = zeros(length(vec_Temperature))
data_cv = zeros(length(vec_Temperature))

for index_Temperature = 1:length(vec_Temperature)
    Temperature = vec_Temperature[index_Temperature]
    print(Temperature, "  ")
    if 2.2 < Temperature < 2.4
        period_sample = 100
    else
        period_sample = 50
    end

    vec_magnetization, vec_energy = stream_tw(Config, Temperature, num_sample, num_sample_burn, period_sample)
    
    data_magnetization[index_Temperature] = mean(abs.(vec_magnetization)) /Length^2
    data_chi[index_Temperature] =  var(abs.(vec_magnetization)) /Length^4
    data_energy[index_Temperature] = mean(vec_energy) /Length^2
    data_cv[index_Temperature] =  var(vec_energy) /Length^4
end    

figure()
scatter(vec_Temperature, data_magnetization, 10, "black")
ylabel("Magnetization")
xlabel("Temperature")
figure()
scatter(vec_Temperature, data_chi, 10, "black")
ylabel("Magnetic susceptibility")
xlabel("Temperature")
figure()
scatter(vec_Temperature, data_energy, 10, "black")
ylabel("Energy")
xlabel("Temperature")
figure()
scatter(vec_Temperature, data_cv, 10, "black")
ylabel("Specific heat")
xlabel("Temperature")

function stream_distance(Config, Temperature, num_sample, num_sample_burn, period_sample)
    Length = size(Config)[1]

    vec_distance = zeros(num_sample)
    
    for time = 1:num_sample_burn
        for iteration = 1:Length^2
            x = mod(rand(Int8),Length) + 1
            y = mod(rand(Int8),Length) + 1

            spin_site = Config[x,y]
            E = 4
            if Config[mod(x,Length)+1,y] == spin_site
                E -= 2
            end
            if Config[mod(x-2,Length)+1,y] == spin_site
                E -= 2
            end
            if Config[x,mod(y,Length)+1] == spin_site
                E -= 2
            end
            if Config[x,mod(y-2,Length)+1] == spin_site
                E -= 2
            end

            if E > 0
                Config[x,y] = !Config[x,y] # Flip
            elseif rand(Float64) < exp(2*E/Temperature)
                Config[x,y] = !Config[x,y] # Flip
            end
        end
    end
    
    Config0 = copy(Config)
    for time = 1:num_sample
        for iteration = 1:Length^2 *period_sample
            x = mod(rand(Int8),Length) + 1
            y = mod(rand(Int8),Length) + 1
            
            spin_site = Config[x,y]
            E = 4
            if Config[mod(x,Length)+1,y] == spin_site
                E -= 2
            end
            if Config[mod(x-2,Length)+1,y] == spin_site
                E -= 2
            end
            if Config[x,mod(y,Length)+1] == spin_site
                E -= 2
            end
            if Config[x,mod(y-2,Length)+1] == spin_site
                E -= 2
            end

            if E > 0
                Config[x,y] = !Config[x,y] # Flip
            elseif rand(Float64) < exp(2*E/Temperature)
                Config[x,y] = !Config[x,y] # Flip
            end
        end

        vec_distance[time] = mean(Config0.!=Config)
    end
    
    return vec_distance
end

function stream_tw_distance(Config, Temperature, num_sample, num_sample_burn, period_sample)
    Length = size(Config)[1]

    vec_distance = zeros(num_sample)
    
    for time = 1:num_sample_burn
        for iteration = 1:Length^2
            for x = 1:Length
                for y = 1:Length
                    spin_site = Config[x,y]
                    E = 4
                    if Config[mod(x,Length)+1,y] == spin_site
                        E -= 2
                    end
                    if Config[mod(x-2,Length)+1,y] == spin_site
                        E -= 2
                    end
                    if Config[x,mod(y,Length)+1] == spin_site
                        E -= 2
                    end
                    if Config[x,mod(y-2,Length)+1] == spin_site
                        E -= 2
                    end

                    if E > 0
                        Config[x,y] = !Config[x,y] # Flip
                    elseif rand(Float64) < exp(2*E/Temperature)
                        Config[x,y] = !Config[x,y] # Flip
                    end
                end
            end
        end
    end
    
    Config0 = copy(Config)
    for time = 1:num_sample
        for iteration = 1:period_sample
            for x = 1:Length
                for y = 1:Length
                    spin_site = Config[x,y]
                    E = 4
                    if Config[mod(x,Length)+1,y] == spin_site
                        E -= 2
                    end
                    if Config[mod(x-2,Length)+1,y] == spin_site
                        E -= 2
                    end
                    if Config[x,mod(y,Length)+1] == spin_site
                        E -= 2
                    end
                    if Config[x,mod(y-2,Length)+1] == spin_site
                        E -= 2
                    end

                    if E > 0
                        Config[x,y] = !Config[x,y] # Flip
                    elseif rand(Float64) < exp(2*E/Temperature)
                        Config[x,y] = !Config[x,y] # Flip
                    end
                end
            end
        end

        vec_distance[time] = mean(Config0.!=Config)
    end
    
    return vec_distance
end

Length = 32
Config = initialize(Length, "random")

Temperature = 2
num_sample = 1000
num_sample_burn = 0
period_sample = 1

for iteration = 1:1
    vec_distance = stream_distance(Config, Temperature, num_sample, num_sample_burn, period_sample)
    plot(vec_distance)
end


Length = 16

Temperature = 2
num_sample = 1000
num_sample_burn = 0
period_sample = 1

for iteration = 1:10
    Config = initialize(Length, "random")
    vec_magnetization, vec_energy = stream_glauber(Config, Temperature, num_sample, num_sample_burn, period_sample)
    plot(vec_magnetization./Length./Length)
end


function autocorrelation(data, k)
    n = length(data)
    mu = mean(data)
    variance = var(data)
    
    corr = 0
    for t = 1:(n - k)
        corr += (data[t] - mu) * (data[t + k] - mu)
    end
    
    return corr / (variance * (n - k))
end
        
function autocorrelation_function(data, kmax)
    corr = zeros(kmax)
    for k = 0:kmax-1
        corr[k+1] = autocorrelation(data, k)
    end
    return corr
end
    
function autocorrelation_time(data, kmax)
    corr = zeros(kmax)
    for k = 0:kmax-1
        corr[k+1] = autocorrelation(data, k)
    end
    if length(findall(corr.<exp(-1))) == 0
        return Inf
    end
    t2 = findall(corr.<exp(-1))[1]
    t1 = t2 -1
    return (t2*(corr[t1]-exp(-1)) + t1*(-corr[t2]+exp(-1))) /(corr[t1]-corr[t2])
end

Length = 32
Config = initialize(Length, "random")

Temperature = 2.3
num_sample = Int(1e5)
num_sample_burn = 2000
period_sample = 1

for iteration = 1:5
    vec_magnetization, vec_energy = stream(Config, Temperature, num_sample, num_sample_burn, period_sample)
    println(autocorrelation_time(abs.(vec_magnetization), 1000))
    semilogy(autocorrelation_function(abs.(vec_magnetization), 1000))
end

ylim([exp(-5), 1])
yticks([1, exp(-1), exp(-2), exp(-3), exp(-4), exp(-5)],["1", L"$1/e$", L"$1/e^2$", L"$1/e^3$", L"$1/e^4$", L"$1/e^5$"])
grid()

Length = 32
Config = initialize(Length, "random")

Temperature = 2.27
num_sample = Int(1e6)
num_sample_burn = 2000
period_sample = 1

vec_magnetization, vec_energy = stream_tw(Config, Temperature, num_sample, num_sample_burn, period_sample)
println(autocorrelation_time(abs.(vec_magnetization), 1000))
semilogy(autocorrelation_function(abs.(vec_magnetization), 1000))

vec_magnetization, vec_energy = stream(Config, Temperature, num_sample, num_sample_burn, period_sample)
println(autocorrelation_time(abs.(vec_magnetization), 1000))
semilogy(autocorrelation_function(abs.(vec_magnetization), 1000))

vec_magnetization, vec_energy = stream_glauber_tw(Config, Temperature, num_sample, num_sample_burn, period_sample)
println(autocorrelation_time(abs.(vec_magnetization), 1000))
semilogy(autocorrelation_function(abs.(vec_magnetization), 1000))

vec_magnetization, vec_energy = stream_glauber(Config, Temperature, num_sample, num_sample_burn, period_sample)
println(autocorrelation_time(abs.(vec_magnetization), 1000))
semilogy(autocorrelation_function(abs.(vec_magnetization), 1000))

yticks([1, exp(-1), exp(-2), exp(-3)],["1", L"$1/e$", L"$1/e^2$", L"$1/e^3$"])
ylim([exp(-3), 1])
grid()

xlabel("Time")
ylabel("Autocorrelation")
legend(["Metropolis, Typewriter", "Metropolis, Ranom-site", "Glauber, Typewriter", "Glauber, Random-site"])

Length = 24
Config = initialize(Length, "random")

Temperature = 2
num_sample = Int(1e5)
num_sample_burn = 100
period_sample = 1

vec_Temperature = vcat(1.5:0.1:2.2, 2.21:0.01:2.49, 2.5:0.1:3)
data_tau_m = zeros(length(vec_Temperature))
data_tau_e = zeros(length(vec_Temperature))
data_tau_m_error = zeros(length(vec_Temperature))
data_tau_e_error = zeros(length(vec_Temperature))

iteration_max = 10
vec_tau_m = zeros(iteration_max)
vec_tau_e = zeros(iteration_max)

for index_Temperature = 1:length(vec_Temperature)
    Temperature = vec_Temperature[index_Temperature]
    for iteration = 1:iteration_max
        vec_magnetization, vec_energy = stream_tw(Config, Temperature, num_sample, num_sample_burn, period_sample)
        vec_tau_m[iteration] = autocorrelation_time(abs.(vec_magnetization), 200)
        vec_tau_e[iteration] = autocorrelation_time(vec_energy, 200)        
    end
    data_tau_m[index_Temperature] += mean(vec_tau_m)
    data_tau_e[index_Temperature] += mean(vec_tau_e)
    data_tau_m_error[index_Temperature] += std(vec_tau_m)/sqrt(iteration_max)
    data_tau_e_error[index_Temperature] += std(vec_tau_e)/sqrt(iteration_max)
    println(Temperature, "  ", data_tau_m[index_Temperature], "  ", data_tau_e[index_Temperature], "  ", data_tau_m_error[index_Temperature], "  ", data_tau_e_error[index_Temperature])
end

writedlm(@sprintf("data/autocorrelation_time_tw_L%d.txt",Length), hcat(vec_Temperature, data_tau_m, data_tau_m_error, data_tau_e, data_tau_e_error))

#errorbar(vec_Temperature, data_tau_m, data_tau_m_error)
#errorbar(vec_Temperature, data_tau_e, data_tau_e_error)
scatter(vec_Temperature, data_tau_m, 5)
scatter(vec_Temperature, data_tau_e, 5)

figure()
for Length=[16 24 32]
    data = readdlm(@sprintf("data/autocorrelation_time_tw_L%d.txt",Length));
    scatter(data[:,1], data[:,2], 10)
end
legend(["L = 16", "L = 24", "L = 32"])
title("Magnetization autocorrelation time")
xlabel("Temperature")

figure()
for Length=[16 24 32]
    data = readdlm(@sprintf("data/autocorrelation_time_tw_L%d.txt",Length));
    scatter(data[:,1], data[:,4], 10)
end
legend(["L = 16", "L = 24", "L = 32"])
title("Energy autocorrelation time")
xlabel("Temperature")

Nu = 1;
Gamma = 1.9;
Tc = 2.27;

figure()
for Length=[16 24 32]
    data = readdlm(@sprintf("data/autocorrelation_time_tw_L%d.txt",Length));
    scatter(((data[:,1].-Tc)./Tc).*Length^(1/Nu), Length^(-Gamma/Nu) .* (data[:,2]), 10)
end
xlim([-5, 10])
legend(["L = 16", "L = 24", "L = 32"])

Nu = 1;
Gamma = 1;
Tc = 2.27;

figure()
for Length=[16 24 32]
    data = readdlm(@sprintf("data/autocorrelation_time_tw_L%d.txt",Length));
    scatter(((data[:,1].-Tc)./Tc).*Length^(1/Nu), Length^(-Gamma/Nu) .* (data[:,4]), 10)
end
xlim([-5, 10])
legend(["L = 16", "L = 24", "L = 32"])
