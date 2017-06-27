#include "result_r.h"
#include "cec_r_utils.h"

SEXP create_R_result(cec_out * out) {
    int m = out->clustering_vector->len;
    int k = out->initial_k;
    int n = out->centers->n;
    int trimmed_size = out->iterations + 1;

    SEXP energy_vector;
    SEXP clusters_number_vector;
    SEXP assignment_vector;
    SEXP covariance_list;
    SEXP centers_matrix;
    SEXP iterations;

    vec_i *cluster_number_trimmed = vec_i_create_from(trimmed_size, out->clusters_number->ar);
    vec_d *energy_trimmed = vec_d_create_from(trimmed_size, out->energy->ar);
    cec_mat *trimmed_centers = cec_matrix_create(k, n);
    for (int i = 0; i < k; i++)
        array_copy(cec_matrix_const_row(out->centers, i), cec_matrix_row(trimmed_centers, i), n);

    PROTECT(energy_vector = allocVector(REALSXP, trimmed_size));
    PROTECT(clusters_number_vector = allocVector(INTSXP, trimmed_size));
    PROTECT(assignment_vector = allocVector(INTSXP, m));
    PROTECT(covariance_list = allocVector(VECSXP, k));
    PROTECT(iterations = allocVector(INTSXP, 1));
    PROTECT(centers_matrix = create_R_matrix(trimmed_centers));

    INTEGER(iterations)[0] = out->iterations;
    vec_d_copy_to(energy_trimmed, REAL(energy_vector));
    vec_i_copy_to(cluster_number_trimmed, INTEGER(clusters_number_vector));
    vec_i_copy_to(out->clustering_vector, INTEGER(assignment_vector));
    for (int i = 0; i < m; i++)
        INTEGER(assignment_vector)[i]++;
    for (int i = 0; i < k; i++) {
        SEXP covariance;
        PROTECT(covariance = create_R_matrix(out->covriances->mats[i]));
        SET_VECTOR_ELT(covariance_list, i, covariance);
    }

    SEXP ret;
    PROTECT(ret = allocList(6));

    SEXP ret_s = ret;

    SETCAR(ret, assignment_vector);
    SET_TAG(ret, install("cluster"));
    ret = CDR(ret);
    SETCAR(ret, centers_matrix);
    SET_TAG(ret, install("centers"));
    ret = CDR(ret);
    SETCAR(ret, energy_vector);
    SET_TAG(ret, install("energy"));
    ret = CDR(ret);
    SETCAR(ret, clusters_number_vector);
    SET_TAG(ret, install("nclusters"));
    ret = CDR(ret);
    SETCAR(ret, covariance_list);
    SET_TAG(ret, install("covariances"));
    ret = CDR(ret);
    SETCAR(ret, iterations);
    SET_TAG(ret, install("iterations"));

    UNPROTECT(k + 7);

    return ret_s;
}
