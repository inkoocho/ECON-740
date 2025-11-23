

using CSV, DataFrames
using Random
using Plots

# path to the CSV located next to this file
csv_path = joinpath(@__DIR__, "ehsw_gpu_results.csv")
@assert isfile(csv_path) "ehsw_gpu_results.csv not found at $csv_path"

# read into a DataFrame and print column names
df = CSV.read(csv_path, DataFrame)
cols = string.(names(df))

########################################################################################
########################################################################################
############# HISTOGRAMS OF THE LAST PERIOD (time == 20000) ############################
########################################################################################
########################################################################################
 
# select last 10000 rows (or fewer if df is smaller)
lastn = min(10000, nrow(df))
sub = df[(nrow(df)-lastn+1):nrow(df), :]

# locate columns (case-insensitive)
price_idx = findfirst(nm -> occursin("price", lowercase(string(nm))), names(df))
pi_idx = findfirst(nm -> lowercase(string(nm)) == "pi", names(df))
# try a looser pi match if exact "pi" not found (but avoid matching "price")
if pi_idx === nothing
    pi_idx = findfirst(nm -> occursin("pi", lowercase(string(nm))) && !occursin("price", lowercase(string(nm))), names(df))
end

@assert price_idx !== nothing "Could not find a 'price' column. Available columns: $(names(df))"
@assert pi_idx !== nothing "Could not find a 'pi' column. Available columns: $(names(df))"

price_col = names(df)[price_idx]
pi_col = names(df)[pi_idx]

prices = sub[!, price_col]
pis = sub[!, pi_col]

# histogram of price
histogram(prices; bins=50,
    title="Histogram of price (last $lastn observations)",
    xlabel=string(price_col), ylabel="Frequency")
savefig(joinpath(@__DIR__, "hist_price_last$(lastn).png"))

# histogram of pi
histogram(pis; bins=50,
    title="Histogram of pi (last $lastn observations)",
    xlabel=string(pi_col), ylabel="Frequency")
savefig(joinpath(@__DIR__, "hist_pi_last$(lastn).png"))


p1 = histogram(prices; bins=50, xlabel=string(price_col), ylabel="Frequency",
    title="Histogram of price (last $lastn observations)", legend=false)
p2 = histogram(pis; bins=50, xlabel=string(pi_col), ylabel="Frequency",
    title="Histogram of pi (last $lastn observations)", legend=false)

combined = plot(p1, p2, layout=(1,2), size=(1200,600))
savefig(joinpath(@__DIR__, "hist_price_pi_last$(lastn).png"))



