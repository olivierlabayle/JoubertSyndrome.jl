module JoubertSyndrome

using DataFrames
using CSV
using DataStructures: SortedDict, Reverse
using StatsBase
using CategoricalArrays
using HypothesisTests
using Plots
using StatsPlots
using FreqTables


DIR_OUT = "/Users/olivierlabayle/Documents/JoubertSyndrome"

DATAPATH = joinpath(DIR_OUT, "AA_these_joubert_descriptif_762021_VF.csv")


PHEN_OF_INTEREST = ("Neonatal breathing dysregulations",
                    "Epilepsy",
                    "Obesity or overweight",
                    "Macrocephaly",
                    "Retinal dystrophy",
                    "Coloboma",
                    "Kidney defect",
                    "Liver defect",
                    "Orofacial anomalies",
                    "Multiple pulmonary infections",
                    "Endocrinal and/or genital features ",
                    "polydactyly",
                    "pure JS")

GENES_OF_INTEREST = ("CC2D2A",
                    "CPLANE1",
                    "AHI1",
                    "CEP290",
                    "TMEM67",
                    "KIAA0586",
                    "INPP5E",
                    "CSPP1",
                    "RPGRIP1L",
                    "OFD1",
                    "MKS1",
                    "NPHP1")


data = DataFrame(CSV.File(DATAPATH; missingstrings=["na"]))


function getfigures(x, y)
    ct = freqtable(x, y)
    odds_phen_if_gen = (ct[2,2]/ct[2,1]) / (ct[1,2]/ct[1,1])
    se = sqrt(sum(1 ./ ct))
    lower_bound = exp(log(odds_phen_if_gen) - 1.96*se)
    upper_bound = exp(log(odds_phen_if_gen) + 1.96*se)
    min_freq = minimum(ct)
    reschi = ChisqTest(ct)
    resfish = FisherExactTest(ct[1,1], ct[1,2], ct[2,1], ct[2,2])
    return min_freq, reschi, resfish, odds_phen_if_gen, lower_bound, upper_bound
end


function gene_phenotype_it(data::DataFrame, gene, phenotype)
    data_restrict = dropmissing(data[:, ["GENE", phenotype]])
    x = categorical(data_restrict.GENE .== gene)
    y = categorical(data_restrict[:, phenotype] .== "oui" )
    return getfigures(x, y)
end


function phenotype_phenotype_it(data::DataFrame, phenotype1, phenotype2)
    data_restrict = dropmissing(data[:, [phenotype1, phenotype2]])
    x = categorical(data_restrict[:, phenotype1] .== "oui")
    y = categorical(data_restrict[:, phenotype2] .== "oui")
    return getfigures(x, y)
end


function gene_phenotype_it_results(data::DataFrame, genes, phenotypes)
    results = DataFrame(Gene = String[], Phenotype = String[], 
                PVAL_CHI2 = Float64[], PVAL_FISH = Float64[], MIN_FREQ = Int[],
                OR = Float64[], OR_LB = Float64[], OR_UB = Float64[], 
                RES_CHI2 = Any[], RES_FISH = Any[])
    for gene in genes
        for phenotype in phenotypes
            min_freq, reschi, resfish, odds_phen_if_gen, lower_bound, upper_bound = gene_phenotype_it(data, gene, phenotype)
            push!(results, 
            (gene, phenotype, pvalue(reschi), pvalue(resfish), min_freq, odds_phen_if_gen, lower_bound, upper_bound, reschi, resfish))
        end
    end
    return results
end

function phenotype_phenotype_it_results(data::DataFrame, phenotypes)
    results = DataFrame(Phenotype1 = String[], Phenotype2 = String[], 
                        PVAL_CHI2 = Float64[], PVAL_FISH = Float64[], 
                        MIN_FREQ = Int[], OR = Float64[], OR_LB = Float64[], OR_UB = Float64[],
                        RES_CHI2 = Any[], RES_FISH = Any[])
    for phenotype1 in phenotypes
        for phenotype2 in phenotypes
            if phenotype1 != phenotype2
                min_freq, reschi, resfish, odds_phen_if_gen, lower_bound, upper_bound = phenotype_phenotype_it(data, phenotype1, phenotype2)
                push!(results, 
                (phenotype1, phenotype2, pvalue(reschi), pvalue(resfish), min_freq, odds_phen_if_gen, lower_bound, upper_bound, reschi, resfish))
            end
        end
    end
    return results
end


gene_phen_results = gene_phenotype_it_results(data, GENES_OF_INTEREST, PHEN_OF_INTEREST)

phen_phen_results = phenotype_phenotype_it_results(data, PHEN_OF_INTEREST)

CSV.write(joinpath(DIR_OUT, "gene_phen.csv"), gene_phen_results[:, Not(["RES_CHI2", "RES_FISH"])])
CSV.write(joinpath(DIR_OUT, "phen_phen.csv"), phen_phen_results[:, Not(["RES_CHI2", "RES_FISH"])])

# genes = categorical(data.GENE)
# plot(countmap(genes), 
#         seriestype=:bar, 
#         xrotation=45, 
#         title="Number of patients for each gene",
#         legend=nothing)



end
