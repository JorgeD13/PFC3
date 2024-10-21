import random
import counter_based_generator
import fortuna_generator
import chaotic_generator
import mersenne_twister
import proposed
import tests

import estimadores_entropia as ee

# n = 100000
# Desde 10000 no hay problema
n = 10000

print("Fortuna")
g2 = fortuna_generator.Accumulator().random_data(n)
print(g2)

print("Chaotic")
g3 = chaotic_generator.mapa_logistico(num_iteraciones=n)
print(g3)

print("mersenne twister")
g4 = mersenne_twister.mersenne_twister(n)
print(g4)

print("Proposed")
g5 = proposed.ProposedAccumulator().random_data(n)
print(g5)

# print("MONOBIT")
# print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
# print(str(tests.monobit(g2)) + "\t\t" + str(tests.monobit(g3)) + "\t\t" + str(tests.monobit(g4)) + "\t\t" + str(tests.monobit(g5)))
#
# print()
#
# print("FREQUENCY BLOCK")
# print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
# print(str(tests.block_frequency(g2)) + "\t\t" + str(tests.block_frequency(g3)) + "\t\t" + str(tests.block_frequency(g4)) + "\t\t" + str(tests.block_frequency(g5)))
#
# print()
#
# print("RUNS")
# print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
# print(str(tests.independent_runs(g2)) + "\t\t" + str(tests.independent_runs(g3)) + "\t\t" + str(tests.independent_runs(g4)) + "\t\t" + str(tests.independent_runs(g5)))
#
# print()
#
# print("RUNS IN A BLOCK")
# print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
# print(str(tests.longest_runs(g2)) + "\t\t" + str(tests.longest_runs(g3)) + "\t\t" + str(tests.longest_runs(g4)) + "\t\t" + str(tests.longest_runs(g5)))
#
# print()
#
# print("MATRIX RANK")
# print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
# print(str(tests.matrix_rank(g2)) + "\t\t" + str(tests.matrix_rank(g3)) + "\t\t" + str(tests.matrix_rank(g4)) + "\t\t" + str(tests.matrix_rank(g5)))
#
# TESTING
print()

print("SPECTRAL")
print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
print(str(tests.spectral(g2)) + "\t\t" + str(tests.spectral(g3)) + "\t\t" + str(tests.spectral(g4)) + "\t\t" + str(tests.spectral(g5)))
#
# print()
#
# print("NON OVERLAPPING")
# print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
# print(str(tests.non_overlapping_patterns(g2)) + "\t\t" + str(tests.non_overlapping_patterns(g3)) + "\t\t" + str(tests.non_overlapping_patterns(g4)) + "\t\t" + str(tests.non_overlapping_patterns(g5)))
#
# print()
#
# print("OVERLAPPING")
# print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
# print(str(tests.overlapping_patterns(g2)) + "\t\t" + str(tests.overlapping_patterns(g3)) + "\t\t" + str(tests.overlapping_patterns(g4)) + "\t\t" + str(tests.overlapping_patterns(g5)))
#
print()

print("UNIVERSAL")
print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
print(str(tests.universal(g2)) + "\t\t" + str(tests.universal(g3)) + "\t\t" + str(tests.universal(g4)) + "\t\t" + str(tests.universal(g5)))

print()
#
# print("LINEAR COMPLEXITY")
# print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
# print(str(tests.linear_complexity(g2)) + "\t\t" + str(tests.linear_complexity(g3)) + "\t\t" + str(tests.linear_complexity(g4)) + "\t\t" + str(tests.linear_complexity(g5)))
#
# print()
#
# print("SERIAL")
# print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
# print(str(tests.serial(g2)) + "\t\t" + str(tests.serial(g3)) + "\t\t" + str(tests.serial(g4)) + "\t\t" + str(tests.serial(g5)))
#
# print()
#
# print("APPROXIMATE ENTROPY")
# print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
# print(str(tests.approximate_entropy(g2)) + "\t\t" + str(tests.approximate_entropy(g3)) + "\t\t" + str(tests.approximate_entropy(g4)) + "\t\t" + str(tests.approximate_entropy(g5)))
#
# Testing
print()

print("CUMULATIVE SUMS")
print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
print(str(tests.cumulative_sums(g2)) + "\t\t" + str(tests.cumulative_sums(g3)) + "\t\t" + str(tests.cumulative_sums(g4)) + "\t\t" + str(tests.cumulative_sums(g5)))
#
# print()
#
# print("RANDOM EXCURSIONS")
# print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
# print(str(tests.random_excursions(g2)) + "\t\t" + str(tests.random_excursions(g3)) + "\t\t" + str(tests.random_excursions(g4)) + "\t\t" + str(tests.random_excursions(g5)))
#
# print()
#
# print("RANDOM EXCURSIONS VARIANT")
# print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
# print(str(tests.random_excursions_variant(g2)) + "\t\t" + str(tests.random_excursions_variant(g3)) + "\t\t" + str(tests.random_excursions_variant(g4)) + "\t\t" + str(tests.random_excursions_variant(g5)))

print()

print("MAXIMUM LIKELIHOOD ENTROPY")
print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
print(str(ee.mle_entropy_binary(g2)) + "\t\t" + str(ee.mle_entropy_binary(g3)) + "\t\t" + str(ee.mle_entropy_binary(g4)) + "\t\t" + str(ee.mle_entropy_binary(g5)))

print()

print("MILLER MADOW ENTROPY")
print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
print(str(ee.miller_madow_entropy_binary(g2)) + "\t\t" + str(ee.miller_madow_entropy_binary(g3)) + "\t\t" + str(ee.miller_madow_entropy_binary(g4)) + "\t\t" + str(ee.miller_madow_entropy_binary(g5)))

print()

print("CHAO SHEN ENTROPY")
print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
print(str(ee.chao_shen_entropy_binary(g2)) + "\t\t" + str(ee.chao_shen_entropy_binary(g3)) + "\t\t" + str(ee.chao_shen_entropy_binary(g4)) + "\t\t" + str(ee.chao_shen_entropy_binary(g5)))

print()

print("BAYESIAN ENTROPY")
print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
print(str(ee.bayesian_entropy_binary(g2)) + "\t\t" + str(ee.bayesian_entropy_binary(g3)) + "\t\t" + str(ee.bayesian_entropy_binary(g4)) + "\t\t" + str(ee.bayesian_entropy_binary(g5)))

print()

print("JAMES STEIN ENTROPY")
print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
print(str(ee.james_stein_entropy_binary(g2)) + "\t\t" + str(ee.james_stein_entropy_binary(g3)) + "\t\t" + str(ee.james_stein_entropy_binary(g4)) + "\t\t" + str(ee.james_stein_entropy_binary(g5)))

print()

print("JACKKNIFE ENTROPY")
print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
print(str(ee.jackknife_entropy_binary(g2)) + "\t\t" + str(ee.jackknife_entropy_binary(g3)) + "\t\t" + str(ee.jackknife_entropy_binary(g4)) + "\t\t" + str(ee.jackknife_entropy_binary(g5)))

print()

print("BEST UPPER BOUND ENTROPY")
print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
print(str(ee.best_upper_bound_entropy_binary(g2)) + "\t\t" + str(ee.best_upper_bound_entropy_binary(g3)) + "\t\t" + str(ee.best_upper_bound_entropy_binary(g4)) + "\t\t" + str(ee.best_upper_bound_entropy_binary(g5)))

print()

print("KERNEL DENSITY ENTROPY")
print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
print(str(ee.kernel_density_entropy_binary(g2)) + "\t\t" + str(ee.kernel_density_entropy_binary(g3)) + "\t\t" + str(ee.kernel_density_entropy_binary(g4)) + "\t\t" + str(ee.kernel_density_entropy_binary(g5)))

print()

print("KNN ENTROPY")
print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
print(str(ee.knn_entropy_binary(g2)) + "\t\t" + str(ee.knn_entropy_binary(g3)) + "\t\t" + str(ee.knn_entropy_binary(g4)) + "\t\t" + str(ee.knn_entropy_binary(g5)))

print()

print("ZHANG ENTROPY")
print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
print(str(ee.zhang_entropy_binary(g2)) + "\t\t" + str(ee.zhang_entropy_binary(g3)) + "\t\t" + str(ee.zhang_entropy_binary(g4)) + "\t\t" + str(ee.zhang_entropy_binary(g5)))

print()

print("NSB ENTROPY")
print("FORTUNA\t\t\t\t\tCHAOTIC\t\t\t\t\tMERSENNE\t\t\t\t\tPROPOSED")
print(str(ee.nsb_entropy_binary(g2)) + "\t\t" + str(ee.nsb_entropy_binary(g3)) + "\t\t" + str(ee.nsb_entropy_binary(g4)) + "\t\t" + str(ee.nsb_entropy_binary(g5)))


