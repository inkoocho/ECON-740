

using CSV, DataFrames
using Random
using Plots

# path to the CSV located next to this file
csv_path = joinpath(@__DIR__, "ehsw_cpu_results.csv")
@assert isfile(csv_path) "ehsw_cpu_results.csv not found at $csv_path"

# read into a DataFrame and print column names
df = CSV.read(csv_path, DataFrame)
cols = string.(names(df))

########################################################################################
########################################################################################
############# HISTOGRAMS OF THE LAST PERIOD (time == 20000) ############################
########################################################################################
########################################################################################


# filter rows with time == 20000
df20000 = filter(row -> row.time == 20000, df)
@assert nrow(df20000) > 0 "No observations with time == 20000"

# separate histograms for price and pi
p1 = histogram(df20000.price; bins=30, xlabel="price", ylabel="Frequency",
               title="Histogram of price (time = 20000)", legend=false)
p2 = histogram(df20000.pi;    bins=30, xlabel="pi",    ylabel="Frequency",
               title="Histogram of pi (time = 20000)",    legend=false)

plot(p1, p2; layout=(1,2), size=(900,400))

########################################################################################
########################################################################################
############# Path of the variables of interest from a random simulation################
########################################################################################
########################################################################################

timecol = :time
vars = [:price, :pi, :beta0, :beta1]



# try to discover a simulation id column (common names) or infer by panel structure
idcol = nothing
tcount = length(unique(df[!, timecol]))
for c in names(df)
    if c == timecol || c in vars
        continue
    end
    if length(unique(df[!, c])) * tcount == nrow(df)
        idcol = c
        break
    end
end

if idcol !== nothing
    sims = unique(df[!, idcol])
    chosen = rand(sims)
    subdf = filter(row -> row[idcol] == chosen, df)
    info = string(idcol, "=", chosen)
else
    # fallback: pick a random combination of the other columns (excluding time and the four vars)
    othercols = setdiff(names(df), [timecol; vars])
    if isempty(othercols)
        subdf = df
        info = "entire dataframe (no id/grouping columns)"
    else
        combos = unique(df[:, othercols])
        rowi = rand(1:nrow(combos))
        combo = combos[rowi, :]
        mask = trues(nrow(df))
        for c in othercols
            mask .= mask .& (df[!, c] .== combo[c])
        end
        subdf = df[mask, :]
        info = "selected unique-combination row $rowi"
    end
end

subdf = sort(subdf, timecol)

p1 = plot(subdf[!, timecol], subdf[!, :price],  xlabel="time", ylabel="price",  label="price",  title="price")
p2 = plot(subdf[!, timecol], subdf[!, :pi],     xlabel="time", ylabel="pi",     label="pi",     title="pi")
p3 = plot(subdf[!, timecol], subdf[!, :beta0],  xlabel="time", ylabel="beta0",  label="beta0",  title="beta0")
p4 = plot(subdf[!, timecol], subdf[!, :beta1],  xlabel="time", ylabel="beta1",  label="beta1",  title="beta1")

plt = plot(p1, p2, p3, p4; layout=(4,1), size=(900,900), link=:x)


outdir = joinpath(@__DIR__, "figures")
mkpath(outdir)

# recreate and save histograms for the last period
h_price = histogram(df20000.price; bins=30, xlabel="price", ylabel="Frequency",
                    title="Histogram of price (time = 20000)", legend=false)
h_pi    = histogram(df20000.pi;    bins=30, xlabel="pi",    ylabel="Frequency",
                    title="Histogram of pi (time = 20000)",    legend=false)
h_combo = plot(h_price, h_pi; layout=(1,2), size=(900,400))

savefig(h_combo, joinpath(outdir, "histograms_last_period.png"))


# save the 4-panel timeseries and each individual timeseries plot
savefig(plt, joinpath(outdir, "random_simulation.png"))

################################################################################
################################################################################
############## SUMMARY STATISTICS FOR THE LAST PERIOD ##########################
################################################################################
################################################################################




using Statistics


pi_vals = collect(skipmissing(df20000.pi))
n_nonmissing = length(pi_vals)
n_missing = count(ismissing, df20000.pi)

if n_nonmissing == 0
    @warn "No non-missing pi values in df20000"
else
    summary_df = DataFrame(
        n = n_nonmissing,
        n_missing = n_missing,
        mean = mean(pi_vals),
        median = median(pi_vals),
        std = std(pi_vals),
        var = var(pi_vals),
        min = minimum(pi_vals),
        q25 = quantile(pi_vals, 0.25),
        q75 = quantile(pi_vals, 0.75),
        max = maximum(pi_vals)
    )

    println("Summary of pi in df20000:")
    show(summary_df, allrows=true, allcols=true)

    # save summary
    savepath = joinpath(outdir, "pi_summary_last_period.csv")
    CSV.write(savepath, summary_df)
end


####################### min=2.75e-29. mean=0.80375. max=1.0 ##########################

