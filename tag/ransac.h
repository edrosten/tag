#ifndef TAG_RANSAC_H_
#define TAG_RANSAC_H_

#include <vector>
#include <algorithm>

namespace tag {

/// @defgroup ransac RANSAC robust estimation
/// This group contains a set of RANSAC implementations to estimate an
/// inlier set from a set of correspondences under a transformation.
/// The functions are both templated on the correspondence data type and
/// the estimator for the transformation.

template <class T>
unsigned int getSample(const std::vector<T>& cdf, const T& p){
    return std::lower_bound(cdf.begin(), cdf.end(), p) - cdf.begin();
}

inline void randomTuple(std::vector<unsigned int>& t, unsigned int bound){
    for (unsigned int i=0; i<t.size(); i++) {
        try_again:
        unsigned int r = rand()%bound;
        for (unsigned int j=0; j<i; j++)
            if (r == t[j])
                goto try_again;
        t[i] = r;
    }
}

template <class T>
void randomTuple(const std::vector<T>& cdf, std::vector<unsigned int>& t, T maxp){
    for (unsigned int i=0; i<t.size(); i++) {
        try_again:
        double x = double(rand()*maxp)/RAND_MAX;
        unsigned int r = std::min((unsigned int)cdf.size()-1, getSample(cdf,(T)x));
        for (unsigned int j=0; j<i; j++)
            if (r == t[j])
                goto try_again;
        t[i] = r;
    }
}

/// basic RANSAC implementation. The function is templated on the correspondence data type
/// and the transformation data type which must conform to the following interface:
/// @code
/// class Estimator {
///     void estimate(const std::vector<Correspondence>& matches);
///     double getSqError(const Correspondence& m);
/// };
/// @endcode
/// see the file @ref ransac_estimators.h for some Estimator classes for various transformations.
/// @param[in] matches a vector of matches of type Correspondence
/// @param[in] dist threshold distance for inliers, the square of dist ist compared with Estimator::getSqError
/// @param[in] sampleSize the number of correspondences for a single hypothesis
/// @param[in] tries the number of hypotheses to test
/// @param[in] bestT the Estimator object used to compute the transformation
/// @param[out] isInlier a vector of bools that is set to describe the inlier set of the winning hypothesis
/// @return the number of inliers for the winning hypothesis
/// @ingroup ransac
template <class Correspondence, class Transform>
int findRANSACInliers(const std::vector<Correspondence>& matches, double dist, unsigned int sampleSize, unsigned int tries, Transform& bestT, std::vector<bool>& isInlier){
    assert(isInlier.size() >= matches.size());
    std::vector<bool> tInlier(isInlier.size());
    Transform T = bestT;
    unsigned int i;
    unsigned int most = 0;
    std::vector<unsigned int> index(sampleSize);
    std::vector<Correspondence> sample(sampleSize);
    for (i=0; i<tries; i++) {
        randomTuple(index, matches.size());
        for (unsigned int j=0; j<sampleSize; j++) {
            sample[j] = matches[index[j]];
        }
        T.estimate(sample);
        unsigned int votes = 0;
        for (unsigned int p=0; p<matches.size() && most+p < matches.size()+votes; p++) {
            double d2 = T.getSqError(matches[p]);
            if (d2 < dist*dist) {
                tInlier[p] = true;
                votes++;
            } else
                tInlier[p] = false;
        }
        if (votes > most) {
            most = votes;
            isInlier = tInlier;
            bestT = T;
            if (most == matches.size())
            break;
        }
    }
    if (most > sampleSize) {
        std::vector<Correspondence> sample(most);
        unsigned int j=0;
        for (i=0;i<matches.size();i++)
            if (isInlier[i])
        sample[j++] = matches[i];
        //bestT.estimate(sample);
    }
    return most;
}

template <class Correspondence, class Transform, class Prob>
int findRANSACInliersGuidedBFS(const std::vector<Correspondence>& matches, const Prob& prob, double dist, unsigned int sampleSize, unsigned int tries, unsigned int perPass, Transform& bestT, std::vector<bool>& isInlier)
{
    assert(isInlier.size() >= matches.size());
    unsigned int i,j;
    std::vector<double> cdf(matches.size());
    double psum = 0;
    for (i=0; i<matches.size(); i++) {
        psum += prob(matches[i]);
        cdf[i] = psum;
    }
    std::vector<Transform> T(tries,bestT);
    std::vector<std::pair<int, unsigned int> > score(tries);
    std::vector<Correspondence> sample(sampleSize);
    for (i=0; i<tries; i++) {
        std::vector<unsigned int> index(sampleSize);
        randomTuple(cdf, index, psum);
        for (j=0; j<sampleSize; j++)
        sample[j] = matches[index[j]];
        T[i].estimate(sample);
        score[i].first = 0;
        score[i].second = i;
    }
    unsigned int curr = 0;
    double distSq = dist*dist;
    int toCut = 1;
    while (curr < matches.size()) {
        if (matches.size() - curr < perPass)
        perPass = matches.size() - curr;
        for (i=0; i<score.size(); i++) {
            const Transform& t = T[score[i].second];
            for (unsigned int p=0; p<perPass; p++) {
                double d2 = t.getSqError(matches[curr+p]);
                if (d2 < distSq)
                    score[i].first++;
            }
        }
        curr += perPass;
        int cutoff = (score.size()*2)/3;//tries-((toCut*tries)>>passes);
        toCut *= 2;
        if (cutoff > 1) {
            nth_element(score.begin(), score.begin()+cutoff, score.end(), std::greater<std::pair<int, unsigned int> >());
            score.resize(cutoff);
        } else
            break;
    }
    bestT = T[max_element(score.begin(), score.end())->second];
    int count = 0;
    for (i=0; i<matches.size(); i++) {
        double d2 = bestT.getSqError(matches[i]);
        if (d2 < distSq) {
            isInlier[i] = true;
            count++;
        } else
            isInlier[i] = false;
    }
    return count;
}

template <class Correspondence, class Transform>
int findRANSACInliersBFS(const std::vector<Correspondence>& matches, double dist, unsigned int sampleSize, unsigned int tries, unsigned int perPass, Transform& bestT, std::vector<bool>& isInlier) {
    assert(isInlier.size() >= matches.size());
    unsigned int i,j;
    std::vector<Transform> T(tries,bestT);
    std::vector<std::pair<int, unsigned int> > score(tries);
    std::vector<Correspondence> sample(sampleSize);
    for (i=0; i<tries; i++) {
        std::vector<unsigned int> index(sampleSize);
        randomTuple(index, matches.size());
        for (j=0; j<sampleSize; j++)
            sample[j] = matches[index[j]];
        T[i].estimate(sample);
        score[i].first = 0;
        score[i].second = i;
    }
    unsigned int curr = 0;
    double distSq = dist*dist;
    int toCut = 1;
    unsigned int passes = tries/perPass + 1;
    while (curr < matches.size()) {
        if (matches.size() - curr < perPass)
            perPass = matches.size() - curr;
        for (i=0; i<score.size(); i++) {
            unsigned int t = score[i].second;
            for (unsigned int p=0; p<perPass; p++) {
                double d2 = T[t].getSqError(matches[curr+p]);
                if (d2 < distSq)
                    score[i].first++;
            }
        }
        curr += perPass;
        int cutoff = tries-((toCut*tries)>>passes);
        toCut *= 2;
        if (cutoff > 1) {
            nth_element(score.begin(), score.begin()+cutoff, score.end(), std::greater<std::pair<int, unsigned int> >());
            score.resize(cutoff);
        } else
            break;
    }
    int bestIndex = max_element(score.begin(), score.end())-score.begin();
    bestT = T[score[bestIndex].second];
    int count = 0;
    for (i=0; i<matches.size(); i++) {
        double d2 = bestT.getSqError(matches[i]);
        if (d2 < distSq) {
            isInlier[i] = true;
            count++;
        } else
            isInlier[i] = false;
    }
    return count;
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
    v.resize(j);
}

} // namespace tag

#endif
