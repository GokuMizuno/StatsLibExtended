namespace internal
{

template <typename T>
double mean(T &x)
{
    using std::begin;
    using std::end;
    mean(begin(x), end(x))
}

template <typename T>
double mean(T &a, const T &b)
{
    return AremeticMean(T &a, const T &b);
}

template <typename T>
double ArimeticMean(T &a, const T &b)
{
    double sum = 0.0;
    //Fastest way to check for zero
    // const uint64_t bits = *(reinterpret_cast<uint64_t*>(&value);
    // const uint32_t uBits = (bits & 0xFFFFFFFF00000000) >> 32;
    // return (uBits + uBits == 0);
    //for (auto i
}


double median(double d_arr[], int size)
{
    double *p_sortedArr = new double[size];
    for (int i=0;i<size;i++)
    {
        p_sortedArr[i] = d_arr[i];
    }
}
// double GetMedian(double daArray[], int iSize) {
//     // Allocate an array of the same size and sort it.
//     double* dpSorted = new double[iSize];
//     for (int i = 0; i < iSize; ++i) {
//         dpSorted[i] = daArray[i];
//     }
//     for (int i = iSize - 1; i > 0; --i) {
//         for (int j = 0; j < i; ++j) {
//             if (dpSorted[j] > dpSorted[j+1]) {
//                 double dTemp = dpSorted[j];
//                 dpSorted[j] = dpSorted[j+1];
//                 dpSorted[j+1] = dTemp;
//             }
//         }
//     }

//     // Middle or average of middle values in the sorted array.
//     double dMedian = 0.0;
//     if ((iSize % 2) == 0) {
//         dMedian = (dpSorted[iSize/2] + dpSorted[(iSize/2) - 1])/2.0;
//     } else {
//         dMedian = dpSorted[iSize/2];
//     }
//     delete [] dpSorted;
//     return dMedian;
// }

// double GetMode(double daArray[], int iSize) {
//     // Allocate an int array of the same size to hold the
//     // repetition count
//     int* ipRepetition = new int[iSize];
//     for (int i = 0; i < iSize; ++i) {
//         ipRepetition[i] = 0;
//         int j = 0;
//         bool bFound = false;
//         while ((j < i) && (daArray[i] != daArray[j])) {
//             if (daArray[i] != daArray[j]) {
//                 ++j;
//             }
//         }
//         ++(ipRepetition[j]);
//     }
//     int iMaxRepeat = 0;
//     for (int i = 1; i < iSize; ++i) {
//         if (ipRepetition[i] > ipRepetition[iMaxRepeat]) {
//             iMaxRepeat = i;
//         }
//     }
//     delete [] ipRepetition;
//     return daArray[iMaxRepeat];
// }

// double GetMean(double daArray[], int iSize) {
//     double dSum = daArray[0];
//     for (int i = 1; i < iSize; ++i) {
//         dSum += daArray[i];
//     }
//     return dSum/iSize;
// }
}