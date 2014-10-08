#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "cec.h"
#include "cov_utils.h"

int cec(struct cec_context * context)
{

    struct cec_matrix * X = context->points;
    struct cec_matrix * C = context->centers;

    /* Init local variables     */

    const int m = X->m;
    const int k = C->m;
    const int n = X->n;
    const int max = context->max_iterations;
    const int min_card = context->min_card;

    energy_function * energy_functions = context->energy_functions;
    struct energy_function_context ** h_contexts =
	    context->energy_function_contexts;

    int _k = k;
    int removed_clusters = 0;
    double energy_sum = 0;

    /* Local arrays and matrices        */
    int * cluster = context->clustering_vector;
    int * clusters_number = context->clusters_number;

    int card[k], removed[k], cluster_map[k];

    double clusters_energy[k];

    double * energy = context->energy;

    struct cec_matrix ** covariance_matrices = context->covriances;

    struct cec_matrix * t_mean_matrix = context->temp_data->t_mean_matrix;
    struct cec_matrix * t_matrix_nn = context->temp_data->t_matrix_nn;
    struct cec_matrix * n_covariance_matrix =
	    context->temp_data->n_covariance_matrix;
    struct cec_matrix ** t_covariance_matrices =
	    context->temp_data->t_covariance_matrices;

    /*
     * For each point find its closest center.
     */
    for (int i = 0; i < m; i++)
    {
	double dist = BIG_DOUBLE;
	for (int j = 0; j < k; j++)
	{
	    double dist_temp = dist2(cec_matrix_row(X, i),
		    cec_matrix_row(C, j), n);
	    if (dist > dist_temp)
	    {
		dist = dist_temp;
		cluster[i] = j;
	    }
	}
    }

    if (max == -1)
    {
	/*
	 * Just return initial cluster assignment
	 */
	context->iterations = -1;
	return NO_ERROR;
    }

    context->iterations = 0;

    /* 
     * Init local arrays
     */
    for (int i = 0; i < k; i++)
    {
	card[i] = 0;
	removed[i] = 0;
	cluster_map[i] = i;
    }

    /*
     *  Set centers matrix to 0.0.
     */
    cec_matrix_set(C, 0.0);

    /* 
     * Set cardinality and center means.
     */
    for (int i = 0; i < m; i++)
    {
	int l = cluster[i];
	card[l]++;
	/*
	 * Sum all data vectors that belong to the same cluster.
	 */
	array_add(cec_matrix_row(C, l), cec_matrix_row(X, i), n);
    }
    for (int i = 0; i < k; i++)
    {
	/*
	 * Check if cardinality meets remove condition.
	 */
	if (card[i] < min_card)
	{
	    removed[i] = 1;
	    array_fill(cec_matrix_row(C, i), NAN, n);
	    removed_clusters++;
	    continue;
	}
	/*
	 * Set the center vector as a mean of all data points in cluster.
	 */
	array_mul(cec_matrix_row(C, i), 1.0 / card[i], n);
    }

    /* 
     * Compute initial covariances of groups.
     */

    /* 
     * Init covariances matrices.
     */
    for (int i = 0; i < k; i++)
    {
	cec_matrix_set(covariance_matrices[i], 0.0);
	cec_matrix_set(t_covariance_matrices[i], 0.0);
    }

    /* Compute covariances      */
    for (int i = 0; i < m; i++)
    {
	int l = cluster[i];
	double t_vec[n];

	/* For each data point compute difference between data point and group mean   */
	array_copy(cec_matrix_row(X, i), t_vec, n);
	array_sub(t_vec, cec_matrix_row(C, l), n);

	/* Sum outer product of difference      */
	cec_vector_outer_product(t_vec, t_matrix_nn, n);
	cec_matrix_add(covariance_matrices[l], t_matrix_nn);
    }

    for (int i = 0; i < k; i++)
    {
	if (removed[i])
	    continue;

	/*
	 * Compute covariances by dividing sum of outer products of differences
	 * between data points and clusters means by cardinality of each group
	 */
	cec_matrix_mul(covariance_matrices[i], 1.0 / card[i]);

	/*
	 * Compute energy of each group
	 */
	double hx = energy_functions[i](h_contexts[i], covariance_matrices[i]);
	if (isnan(hx))
	    return *(h_contexts[i]->last_error);

	clusters_energy[i] = compute_energy(m, hx, card[i]);

	/*
	 * Sum of energy of all groups
	 */
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

    /* Set number of clusters at "iteration zero".      */
    clusters_number[0] = _k;

    energy[0] = energy_sum;

    /*
     * If cluster was removed before first iteration - we must handle it.
     */
    int handle_removed_flag = (k == _k) ? 0 : 1;

    /*
     * Main loop.
     */

    for (int iter = (handle_removed_flag ? -1 : 0); iter < max; iter++)
    {
	/*
	 *  Flag that indicates transfer of point.
	 */
	int transfer_flag = 0;
	int removed_last_iteration_flag = 0;

	/* 
	 * Iterate over data points.
	 */
	for (int i = 0; i < m; i++)
	{

	    /*
	     * Group that the point is assigned to.
	     */
	    int l = cluster[i];

	    if (handle_removed_flag && !removed[l])
		continue;

	    /*
	     * Energy of group l after remove point i.
	     */
	    double n_l_energy = NAN;

	    /*
	     * Energy gain.
	     */
	    double energy_gain;

	    /*
	     * New mean of group l after remove point i.
	     */
	    double n_mean[n];

	    if (!removed[l])
	    {

		/*
		 * Compute mean of group 'l' after removing data point 'i' and store
		 * its value in the n_mean.
		 */
		mean_remove_point(cec_matrix_row(C, l), n_mean,
			cec_matrix_row(X, i), card[l], n);

		/*
		 * Compute covariance of group 'l' after removing data point 'i' and store
		 * its value in n_covariance_matrix.
		 */
		cec_cov_remove_point(covariance_matrices[l],
			n_covariance_matrix, cec_matrix_row(C, l),
			cec_matrix_row(X, i), card[l], t_matrix_nn);

		/*
		 * Compute energy of group 'l' after removing data point 'i'.
		 */
		double n_l_hx = energy_functions[l](h_contexts[l], n_covariance_matrix);
		if (isnan(n_l_hx))
		    return *(h_contexts[l]->last_error);

		n_l_energy = compute_energy(m, n_l_hx, card[l] - 1);

		energy_gain = 0;

	    } else
	    {
		energy_gain = INFINITY;
	    }

	    /*
	     * Index of the best group (energy gain) for point 'i' to transfered.
	     */
	    int idx = -1;

	    /*
	     * Energy of the best group after adding point 'i'.
	     */
	    double n_energy;

	    /*
	     * Iterate over clusters that are not removed.
	     */
	    for (int _j = 0; _j < _k; _j++)
	    {

		/*
		 * Get a cluster number from mapping array.
		 */
		int j = cluster_map[_j];

		if ((j == l) || (removed[j] == 1))
		    continue;

		/*
		 * Compute mean of group 'j' after adding data point 'i' and store
		 * its value in the matrix of temporary means t_mean_matrix at 'j'-th row.
		 */
		mean_add_point(cec_matrix_row(C, j),
			cec_matrix_row(t_mean_matrix, j), cec_matrix_row(X, i),
			card[j], n);

		/*
		 * Compute covariance of group 'j' after adding data point 'i' and store
		 * its value in the array of temporary covariance matrices at index 'j'.
		 */
		cec_cov_add_point(covariance_matrices[j],
			t_covariance_matrices[j], cec_matrix_row(C, j),
			cec_matrix_row(X, i), card[j], t_matrix_nn);

		/*
		 * Compute energy of group 'j' after adding data point 'i'.
		 */
		double t_hx = energy_functions[j](h_contexts[j], t_covariance_matrices[j]);
		if (isnan(t_hx))
		    return *(h_contexts[j]->last_error);

		double t_energy = compute_energy(m, t_hx, card[j] + 1);

		if (removed[l] == 1)
		{
		    /*
		     * Since the energy of cluster 'l' was subtracted from energy sum (when 'l' was removed),
		     * gain is only a change of the energy of the cluster 'j' by adding point 'i' to it.
		     */
		    double gain = (t_energy - clusters_energy[j]);
		    if (gain < energy_gain)
		    {
			idx = j;
			energy_gain = gain;
			n_energy = t_energy;
		    }
		} else
		{
		    double gain = (n_l_energy + t_energy)
			    - (clusters_energy[l] + clusters_energy[j]);
		    if (gain < energy_gain)
		    {
			idx = j;
			energy_gain = gain;
			n_energy = t_energy;
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
		clusters_energy[idx] = n_energy;
		cec_matrix_copy_data(t_covariance_matrices[idx],
			covariance_matrices[idx]);
		array_copy(cec_matrix_row(t_mean_matrix, idx),
			cec_matrix_row(C, idx), n);

		if (!removed[l])
		{
		    cec_matrix_copy_data(n_covariance_matrix,
			    covariance_matrices[l]);
		    clusters_energy[l] = n_l_energy;
		    card[l]--;
		    array_copy(n_mean, cec_matrix_row(C, l), n);

		    /*
		     * Check if cluster 'l' needs to be removed.
		     */
		    if (card[l] < min_card)
		    {
			removed_last_iteration_flag = 1;
			/*
			 * If cluster is being removed, subtract its energy from energy sum.
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
	     *  Change mapping and k for removed clusters.
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
	context->iterations = iter + 1;

	/*
	 * If there was no transfer, finish execution.
	 */
	if (!transfer_flag)
	    return NO_ERROR;

	/*
	 * If any cluster is removed in this iteration, we need to perform another iteration
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

