#include "r_result.h"
#include "r_utils.h"

using cec::r::put;

SEXP cec::create_R_result(const clustering_results &out) {
    int m = out.assignment.size();
    int k = out.centers.m;

    SEXP energy_vector;
    SEXP clusters_number_vector;
    SEXP assignment_vector;
    SEXP covariance_list;
    SEXP centers_matrix;
    SEXP iterations;

    PROTECT(energy_vector = put(out.energy));
    PROTECT(clusters_number_vector = put(out.cluster_number));
    PROTECT(assignment_vector = put(out.assignment));
    PROTECT(covariance_list = allocVector(VECSXP, k));
    PROTECT(iterations = put(out.iterations));
    PROTECT(centers_matrix = put(out.centers));

   int *assignment_vector_data = INTEGER(assignment_vector);
   for (int i = 0; i < m; i++) {
       assignment_vector_data[i]++;
   }
    for (int i = 0; i < k; i++) {
        SEXP covariance;
        PROTECT(covariance = put(out.covariances[i]));
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
