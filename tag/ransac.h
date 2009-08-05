#ifndef TAG_RANSAC_H_
#define TAG_RANSAC_H_

#include <vector>
#include <algorithm>

namespace tag {

#ifdef WIN32
///	taken from http://www.azillionmonkeys.com/qed/random.html
inline double drand48(void) {
	const double RS_SCALE  = (1.0 / (1.0 + RAND_MAX));

    double d;
	do {
		d = (((rand () * RS_SCALE) + rand ()) * RS_SCALE + rand ()) * RS_SCALE;
    } while (d >= 1); /* Round off */
    return d;
}
#endif

/// @defgroup ransac RANSAC robust estimation
/// This group contains a set of RANSAC implementations to estimate an
/// inlier set from a set of correspondences under a transformation.
/// The functions are both templated on the correspondence data type and
/// the estimator for the transformation.

 template <class T, class I> void randomTuple(const std::vector<T>& cdf, std::vector<I>& t, T maxp) {
     for (size_t i=0; i<t.size(); i++) {
     try_again:
	 double x = drand48()* maxp;
	 size_t r = std::min(cdf.size()-1, (size_t)(std::lower_bound(cdf.begin(), cdf.end(), (T)x) - cdf.begin()));
	 for (size_t j=0; j<i; j++)
	     if (r == t[j])
		 goto try_again;
	 t[i] = r;
     }
 }

 template <class T> void randomTuple(T& t, unsigned int bound)
 {
     for (size_t i=0; i<t.size(); i++) {
     try_again:
	 size_t r = (size_t)(drand48() * bound);
	 for (size_t j=0; j<i; j++)
	     if (r == t[j])
		 goto try_again;
	 t[i] = r;
     }
 }


/// basic RANSAC implementation. The function is templated on the observation data type
/// and the transformation data type which must conform to the following interface:
/// @code
/// class Estimator {
///     Estimator();
///     // Estimate from a sequence of observations
///     template <class It> bool estimate(It begin, It End);
///     // Check whether the given observation is an inlier for this estimate (with specified tolerance)
///     template <class Obs, class Tol> bool isInlier(const Obs& obs, const Tol& tolerance) const;
///     // the number of observations to estimate one hypothesis
///     static const int hypothesis_size = XXX;
/// };
/// @endcode
/// see the file @ref ransac_estimators.h for some Estimator classes for various transformations.
/// @param[in] observations a vector of observations (usually point matches)
/// @param[in] sample_size the number of samples used estimate a transformation
/// @param[in] tolerance the tolerance (passed with each observation) to the transformation to check for inliers
/// @param[in] N the number of hypotheses to test
/// @param[out] best the transformation hypothesis with the highest inlier count
/// @param[out] inlier a vector of bools that describes the inlier set of the winning hypothesis
/// @return the number of inliers for the winning hypothesis
/// @ingroup ransac

template <class Obs, class Trans, class Tol> size_t find_RANSAC_inliers(const std::vector<Obs>& observations, const Tol& tolerance, size_t N,
									Trans& best, std::vector<bool>& inlier, int sample_size = Trans::hypothesis_size )
{
    std::vector<bool> thisInlier(observations.size());
    size_t bestScore = 0;
    std::vector<size_t> sample_index(sample_size);
    std::vector<Obs> sample(sample_size);
    while (N--) {
	randomTuple(sample_index, observations.size());
	for (int i=0;i<sample_size; i++)
	    sample[i] = observations[sample_index[i]];
	Trans thisT(best);
	thisT.estimate(sample.begin(), sample.end());
	size_t score = 0;
	for (size_t i=0; i<observations.size(); i++) {
	    const Obs& o = observations[i];
	    if (thisT.isInlier(o, tolerance)) {
		thisInlier[i] = true;
		score++;
	    } else
		thisInlier[i] = false;
	}
	if (score > bestScore) {
	    bestScore = score;
	    inlier = thisInlier;
	    best = thisT;
	}
    }
    return bestScore;
}

/// backwards compatibility interface to find_RANSAC_inliers
/// @deprecated
/// @ingroup ransac
template <class Obs, class Trans, class Tol> size_t find_RANSAC_inliers(const std::vector<Obs>& observations, int sample_size, const Tol& tolerance, size_t N,
									Trans& best, std::vector<bool>& inlier)
{
    return find_RANSAC_inliers(observations, tolerance, N, best, inlier, sample_size);
}

/// basic MSAC implementation. The function is templated on the observation data type
/// and the transformation data type which must conform to the following interface:
/// @code
/// class Estimator {
///     Estimator();
///     // Estimate from a sequence of observations
///     template <class It> bool estimate(It begin, It End);
///     // return the score for the given observation
///     template <class Obs, class Tol> double score(const Obs& obs) const;
///     // the number of observations to estimate one hypothesis
///     static const int hypothesis_size = XXX;
/// };
/// @endcode
/// see the file @ref ransac_estimators.h for some Estimator classes for various transformations.
/// @param[in] observations a vector of observations (usually point matches)
/// @param[in] sample_size the number of samples used estimate a transformation
/// @param[in] tolerance the tolerance (passed with each observation) to the transformation to check for inliers
/// @param[in] N the number of hypotheses to test
/// @param[out] best the transformation hypothesis with the highest inlier count
/// @param[out] inlier a vector of bools that describes the inlier set of the winning hypothesis
/// @return the score of the winning hypothesis
/// @ingroup ransac
template <class Obs, class Trans, class Tol> double find_MSAC_inliers(const std::vector<Obs>& observations, const Tol& tolerance, size_t N,
									Trans& best, std::vector<bool>& inlier, int sample_size = Trans::hypothesis_size)
{
    std::vector<bool> thisInlier(observations.size());
    const double toleranceSquared = tolerance * tolerance;
    double bestScore = observations.size() * toleranceSquared;

    std::vector<size_t> sample_index(sample_size);
    std::vector<Obs> sample(sample_size);
    while (N--) {
	randomTuple(sample_index, observations.size());
	for (int i=0;i<sample_size; i++)
	    sample[i] = observations[sample_index[i]];
	Trans thisT(best);
	thisT.estimate(sample.begin(), sample.end());
	double score = 0;
	for (size_t i=0; i<observations.size(); i++) {
	    const Obs& o = observations[i];
	    const double s = thisT.score(o);
	    if (s < toleranceSquared) {
		thisInlier[i] = true;
		score += s;
	    } else {
		thisInlier[i] = false;
		score += toleranceSquared;
	    }
	}
	if (score < bestScore) {
	    bestScore = score;
	    inlier = thisInlier;
	    best = thisT;
	}
    }
    return bestScore;
}

/// backwards compatibility interface to find_MSAC_inliers
/// @deprecated
/// @ingroup ransac
template <class Obs, class Trans, class Tol> double find_MSAC_inliers(const std::vector<Obs>& observations, int sample_size, const Tol& tolerance, size_t N,
									Trans& best, std::vector<bool>& inlier)
{
    return find_MSAC_inliers(observations, tolerance, N, best, inlier, sample_size);
}


inline double getShrinkRatio(unsigned int H, unsigned int N, unsigned int B)
{
     return pow(double(H), -double(B)/N);
}

/// Guided breadth-first RANSAC implementation. The function is templated on the observation data type,
/// the probability density for observation correctness, the tolerance for inliers,
/// and the transformation data type which must conform to the following interface:
/// @code
/// class Estimator {
///     Estimator();
///     // Estimate from a sequence of observations
///     template <class It> bool estimate(It begin, It End);
///     // Check whether the given observation is an inlier for this estimate (with specified tolerance)
///     template <class Obs, class Tol> bool isInlier(const Obs& obs, const Tol& tolerance) const;
///     // the number of observations to estimate one hypothesis
///     static const int hypothesis_size = XXX;
/// };
/// @endcode
/// All hypotheses are generated first, and preemptively discarded as more observations are examined.
/// The sample sets for hypotheses are drawn probabilistically using the specified probability density.
/// see the file @ref ransac_estimators.h for some Estimator classes for various transformations.
/// @param[in] observations a vector of observations (usually point matches)
/// @param[in] prob the probability predicate applied (as a function object) to each observation to compute its prior probability of being correct
/// @param[in] sample_size the number of samples used estimate a transformation
/// @param[in] tolerance the tolerance (passed with each observation) to the transformation to check for inliers
/// @param[in] N the number of hypotheses to test
/// @param[in] block_size the number of hypotheses to test in a block (between culling hypotheses)
/// @param[out] best the transformation hypothesis with the highest inlier count
/// @param[out] inlier a vector of bools that describes the inlier set of the winning hypothesis
/// @return the number of inliers for the winning hypothesis
/// @ingroup ransac
template <class Obs, class Trans, class Tol, class Prob>
size_t find_RANSAC_inliers_guided_breadth_first(const std::vector<Obs>& observations, const Prob& prob, const Tol& tolerance, size_t N, size_t block_size,
							Trans& best, std::vector<bool>& inlier, int sample_size = Trans::hypothesis_size)
    {
	std::vector<Trans> hypotheses(N,best);
	std::vector<std::pair<int,size_t> > score(N);
	std::vector<size_t> sample_index(sample_size);
	std::vector<Obs> sample(sample_size);
	std::vector<double> cdf(observations.size());
	cdf[0] = prob(observations[0]);
	for (size_t i=1; i<observations.size(); ++i)
	    cdf[i] = cdf[i-1] + prob(observations[i]);
	const double psum = cdf.back();

	for (size_t i=0; i<hypotheses.size(); i++) {
	    do {
		randomTuple(cdf, sample_index, psum);
		for (int s=0; s<sample_size; ++s)
		    sample[s] = observations[sample_index[s]];
	    }
	    while (!hypotheses[i].estimate(sample.begin(), sample.end()));
	    score[i] = std::make_pair(0,i);
	}
	size_t m = 0;
	const double factor = getShrinkRatio(N, observations.size(), block_size);
	while (m < observations.size()) {
	    size_t end = std::min(observations.size(), m+block_size);
	    for (size_t i=0; i<score.size(); i++) {
		const Trans& thisT = hypotheses[score[i].second];
		size_t s = 0;
		for (size_t j=m; j!=end; j++) {
		    if (thisT.isInlier(observations[j], tolerance))
			++s;
		}
		score[i].first += s;
	    }
	    unsigned int cutoff = (unsigned int)(score.size() * factor);
	    if (cutoff == 0)
		break;
	    std::nth_element(score.begin(), score.end(), score.begin()+cutoff, std::greater<std::pair<int,size_t> >());
	    score.resize(cutoff);
	    m = end;
	}
	size_t best_index = std::max_element(score.begin(), score.end())->second;
	best = hypotheses[best_index];
	size_t count = 0;
	inlier.resize(observations.size());
	for (size_t i=0; i<observations.size(); i++) {
	    if (best.isInlier(observations[i], tolerance)) {
		inlier[i] = true;
		++count;
	    } else
		inlier[i] = false;
	}
	return count;
    }

/// backwards compatibility interface to find_RANSAC_inliers_guided_breadth_first
/// @deprecated
/// @ingroup ransac
template <class Obs, class Trans, class Tol, class Prob>
size_t find_RANSAC_inliers_guided_breadth_first(const std::vector<Obs>& observations, const Prob& prob, int sample_size, const Tol& tolerance, size_t N, size_t block_size,
							Trans& best, std::vector<bool>& inlier)
{
    return find_RANSAC_inliers_guided_breadth_first(observations, prob, tolerance, N, block_size, best, inlier, sample_size);
}

/// Breadth-first RANSAC implementation. The function is templated on the observation data type, the tolerance for inliers,
/// and the transformation data type which must conform to the following interface:
/// @code
/// class Estimator {
///     Estimator();
///     // Estimate from a sequence of observations
///     template <class It> bool estimate(It begin, It End);
///     // Check whether the given observation is an inlier for this estimate (with specified tolerance)
///     template <class Obs, class Tol> bool isInlier(const Obs& obs, const Tol& tolerance) const;
///     // the number of observations to estimate one hypothesis
///     static const int hypothesis_size = XXX;
/// };
/// @endcode
/// All hypotheses are generated first, and preemptively discarded as more observations are examined.
/// see the file @ref ransac_estimators.h for some Estimator classes for various transformations.
/// @param[in] observations a vector of observations (usually point matches)
/// @param[in] sample_size the number of samples used estimate a transformation
/// @param[in] tolerance the tolerance (passed with each observation) to the transformation to check for inliers
/// @param[in] N the number of hypotheses to test
/// @param[in] block_size the number of hypotheses to test in a block (between culling hypotheses)
/// @param[out] best the transformation hypothesis with the highest inlier count
/// @param[out] inlier a vector of bools that describes the inlier set of the winning hypothesis
/// @return the number of inliers for the winning hypothesis
/// @ingroup ransac
template <class Obs, class Trans, class Tol> size_t find_RANSAC_inliers_breadth_first(const std::vector<Obs>& observations, const Tol& tolerance, size_t N, size_t block_size,
											  Trans& best, std::vector<bool>& inlier, int sample_size = Trans::hypothesis_size)
    {
	std::vector<Trans> hypotheses(N,best);
	std::vector<std::pair<int,size_t> > score(N);
	std::vector<size_t> sample_index(sample_size);
	std::vector<Obs> sample(sample_size);
	for (size_t i=0; i<hypotheses.size(); i++) {
	    do {
		randomTuple(sample_index, observations.size());
		for (int s=0; s<sample_size; ++s)
		    sample[s] = observations[sample_index[s]];
	    }
	    while (!hypotheses[i].estimate(sample.begin(), sample.end()));
	    score[i] = std::make_pair(0,i);
	}
	size_t m = 0;
	const double factor = getShrinkRatio(N, observations.size(), block_size);
	while (m < observations.size()) {
	    size_t end = std::min(observations.size(), m+block_size);
	    for (size_t i=0; i<score.size(); i++) {
		const Trans& thisT = hypotheses[score[i].second];
		size_t s = 0;
		for (size_t j=m; j!=end; j++) {
		    if (thisT.isInlier(observations[j], tolerance))
			++s;
		}
		score[i].first += s;
	    }
	    unsigned int cutoff = (unsigned int)(score.size() * factor);
	    if (cutoff == 0)
		break;
	    std::nth_element(score.begin(), score.end(), score.begin()+cutoff, std::greater<std::pair<int,size_t> >());
	    score.resize(cutoff);
	    m = end;
	}
	size_t best_index = std::max_element(score.begin(), score.end())->second;
	best = hypotheses[best_index];
	size_t count = 0;
	inlier.resize(observations.size());
	for (size_t i=0; i<observations.size(); i++) {
	    if (best.isInlier(observations[i], tolerance)) {
		inlier[i] = true;
		++count;
	    } else
		inlier[i] = false;
	}
	return count;
    }

/// backwards compatibility interface to find_RANSAC_inliers_breadth_first
/// @deprecated
/// @ingroup ransac
template <class Obs, class Trans, class Tol> size_t find_RANSAC_inliers_breadth_first(const std::vector<Obs>& observations, int sample_size, const Tol& tolerance, size_t N, size_t block_size,
											  Trans& best, std::vector<bool>& inlier)
{
    return find_RANSAC_inliers_breadth_first(observations, tolerance, N, block_size, best, inlier, sample_size);
}

/// helper function for use with ransac functions. It removes elements
/// of an input vector based on another vector which is interpreted as
/// boolean flags signaling removal from the first vector or not.
/// @param[in,out] v vector with data to remove
/// @param[in] flag vector interpreted as flags, if false the corresponding element from v is removed
/// @ingroup ransac
template <class T, class U>
void remove_if_false(std::vector<T>& v, const std::vector<U>& flag) {
    assert(v.size() <= flag.size());
    unsigned int i=0,j;
    while (i<v.size() && flag[i])
        i++;
    j = i++;
    for (;i<v.size();i++) {
        if (flag[i])
            v[j++] = v[i];
    }
    v.erase(v.begin() + j, v.end());
}

} // namespace tag

#endif
