#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "cec.h"
#include "cov_utils.h"
#include "model.h"

int cec(struct cec_context * context)
{
    const struct cec_matrix * X = context->input->points;

    const int m = X->m;
    const int k = context->input->centers->m;
    const int n = X->n;
    const int max = context->input->max_iterations;
    const int min_card = context->input->min_card;

    struct cec_model ** models = context->input->models;

    int * cluster = context->results->clustering_vector;
    int * clusters_number = context->results->clusters_number;
    struct cec_matrix ** covariance_matrices = context->results->covriances;
    double * energy = context->results->energy;
    
    struct cec_matrix * C = context->results->centers;
    cec_matrix_copy_data(context->input->centers, C);

    struct cec_matrix * t_mean_matrix = context->temp_data->t_mean_matrix;
    
    int _k = k;
    int removed_clusters = 0;
    double energy_sum = 0;
    int card[k], removed[k], cluster_map[k];
    double clusters_energy[k];

    struct cec_matrix * t_matrix_nn = context->temp_data->t_matrix_nn;
    struct cec_matrix * n_covariance_matrix =
            context->temp_data->n_covariance_matrix;
    struct cec_matrix ** t_covariance_matrices =
            context->temp_data->t_covariance_matrices;

    context->results->iterations = 0;

    /*
     * Assign points to its closest clusters and calculate clusters means.
     */

    for (int i = 0; i < m; i++)
    {
        double dist = BIG_DOUBLE;
        for (int j = 0; j < k; j++)
        {
            double dist_temp = dist2(cec_matrix_const_row(X, i),
                    cec_matrix_const_row(C, j), n);
            if (dist > dist_temp)
            {
                dist = dist_temp;
                cluster[i] = j;
            }
        }
    }

    for (int i = 0; i < k; i++)
    {
        card[i] = 0;
        removed[i] = 0;
        cluster_map[i] = i;
    }

    cec_matrix_set(C, 0.0);

    for (int i = 0; i < m; i++)
    {
        int l = cluster[i];
        card[l]++;
        array_add(cec_matrix_row(C, l), cec_matrix_const_row(X, i), n);
    }

    for (int i = 0; i < k; i++)
    {
        if (card[i] < min_card)
        {
            removed[i] = 1;
            array_fill(cec_matrix_row(C, i), NAN, n);
            removed_clusters++;
            continue;
        }
        array_mul(cec_matrix_row(C, i), 1.0 / card[i], n);
    }

    /* 
     * Compute initial covariances using maximum likelihood estimator.
     */

    for (int i = 0; i < k; i++)
    {
        cec_matrix_set(covariance_matrices[i], 0.0);
        cec_matrix_set(t_covariance_matrices[i], 0.0);
    }

    for (int i = 0; i < m; i++)
    {
        int l = cluster[i];
        double t_vec[n];

        array_copy(cec_matrix_const_row(X, i), t_vec, n);
        array_sub(t_vec, cec_matrix_const_row(C, l), n);

        cec_vector_outer_product(t_vec, t_matrix_nn, n);
        cec_matrix_add(covariance_matrices[l], t_matrix_nn);
    }

    for (int i = 0; i < k; i++)
    {
        if (removed[i])
            continue;

        cec_matrix_mul(covariance_matrices[i], 1.0 / card[i]);

        double energy = cluster_energy(models[i], covariance_matrices[i], card[i], m);
        
        if (isnan(energy))
            return cluster_energy_get_last_error(models[i]);

        clusters_energy[i] = energy;

        energy_sum += clusters_energy[i];
    }

    /* 
     * Change cluster mapping for removed clusters.
     */
    for (int i = k; i > 0; i--)
    {
        if (removed[i - 1])
        {
            cluster_map[i - 1] = cluster_map[_k - 1];
            _k--;
        }
    }

    if (_k == 0)
    {
        return ALL_CLUSTERS_REMOVED_ERROR;
    }

    clusters_number[0] = _k;

    energy[0] = energy_sum;

    /*
     * Special case when a cluster was removed before the first iteration.
     */
    int handle_removed_flag = (k == _k) ? 0 : 1;

    /*
     * Main loop.
     */
    for (int iter = (handle_removed_flag ? -1 : 0); iter < max; iter++)
    {
        int transfer_flag = 0;
        int removed_last_iteration_flag = 0;

        for (int i = 0; i < m; i++)
        {

            int l = cluster[i];

            if (handle_removed_flag && !removed[l])
                continue;

            /*
             * Energy and mean vector of cluster 'l' after removing point 'i'.
             * Initializing to NAN to get rid of compiler warning.
             */
            double n_l_energy = NAN;
            double n_l_mean[n];

            double energy_gain;

            if (!removed[l])
            {
                /*
                 * Compute mean of group 'l' after removing data point 'i' and store
                 * its value in the n_mean.
                 */
                mean_remove_point(cec_matrix_const_row(C, l), n_l_mean, cec_matrix_const_row(X, i), card[l], n);

                /*
                 * Compute covariance of group 'l' after removing data point 'i' and store
                 * its value in n_covariance_matrix.
                 */
                cec_cov_remove_point(covariance_matrices[l], n_covariance_matrix, cec_matrix_const_row(C, l),
                        cec_matrix_const_row(X, i), card[l], t_matrix_nn);

                /*
                 * Compute energy of group 'l' after removing data point 'i'.
                 */
                
                double energy = cluster_energy(models[l], n_covariance_matrix, card[l] - 1, m);
        
                if (isnan(energy))
                    return cluster_energy_get_last_error(models[l]);
                
                n_l_energy = energy;

                energy_gain = 0;
            } 
            else
            {
                energy_gain = INFINITY;
            }

            int idx = -1;

            double best_energy;

            for (int _j = 0; _j < _k; _j++)
            {
                int j = cluster_map[_j];

                if ((j == l) || (removed[j] == 1))
                    continue;

                mean_add_point(cec_matrix_row(C, j), cec_matrix_row(t_mean_matrix, j), cec_matrix_const_row(X, i),
                        card[j], n);

                cec_cov_add_point(covariance_matrices[j], t_covariance_matrices[j], cec_matrix_const_row(C, j),
                        cec_matrix_const_row(X, i), card[j], t_matrix_nn);

                double energy = cluster_energy(models[j], t_covariance_matrices[j], card[j] + 1, m);
        
                if (isnan(energy))
                    return cluster_energy_get_last_error(models[j]);
                 
                double t_energy = energy;

                if (removed[l] == 1)
                {
                    /*
                     * Since the energy of cluster 'l' was subtracted from the energy sum (when 'l' was removed),
                     * gain is only the change of the energy of cluster 'j' by adding point 'i'.
                     */
                    double gain = (t_energy - clusters_energy[j]);
                    if (gain < energy_gain)
                    {
                        idx = j;
                        energy_gain = gain;
                        best_energy = t_energy;
                    }
                } else
                {
                    double gain = (n_l_energy + t_energy)
                            - (clusters_energy[l] + clusters_energy[j]);
                    if (gain < energy_gain)
                    {
                        idx = j;
                        energy_gain = gain;
                        best_energy = t_energy;
                    }
                }
            }

            /* 
             * Transfer point 'i' to cluster 'idx'.
             */
            if (idx != -1)
            {
                cluster[i] = idx;
                card[idx]++;
                clusters_energy[idx] = best_energy;
                cec_matrix_copy_data(t_covariance_matrices[idx], covariance_matrices[idx]);
                array_copy(cec_matrix_const_row(t_mean_matrix, idx), cec_matrix_row(C, idx), n);

                if (!removed[l])
                {
                    cec_matrix_copy_data(n_covariance_matrix, covariance_matrices[l]);
                    clusters_energy[l] = n_l_energy;
                    card[l]--;
                    array_copy(n_l_mean, cec_matrix_row(C, l), n);

                    if (card[l] < min_card)
                    {
                        removed_last_iteration_flag = 1;

                        /*
                         * If the cluster is being removed, subtract its energy from the energy sum.
                         */
                        energy_sum -= clusters_energy[l];
                        array_fill(cec_matrix_row(C, l), NAN, n);
                        cec_matrix_set(covariance_matrices[l], NAN);
                        removed[l] = 1;
                    }
                }
                energy_sum += energy_gain;
                transfer_flag = 1;
            }

            /*
             *  Change the mapping and _k for removed clusters.
             */
            for (int _j = _k; _j > 0; _j--)
            {
                int j = cluster_map[_j - 1];
                if (removed[j])
                {
                    cluster_map[_j - 1] = cluster_map[_k - 1];
                    _k--;
                }
            }
        }

        energy[iter + 1] = energy_sum;
        clusters_number[iter + 1] = _k;
        context->results->iterations = iter + 1;

        if (!transfer_flag)
            return NO_ERROR;

        /*
         * If cluster was removed in this iteration, we need to perform another iteration
         * considering points that are not assigned. It will prevent energy dips.
         */

        if (removed_last_iteration_flag)
        {
            handle_removed_flag = 1;
            iter--;
        } else
            handle_removed_flag = 0;
    }
    return NO_ERROR;
}

