/*G. Stangler
  Here we find the mean
  We compute and return the arithemetic mean, the geometric mean, and the harmonic mean
*/

/*template <class FWIT>
void func(FWIT a, const FWIT b)
{
    while (a != b)
    {
        std::cout << "Val:" << *a << std::endl;
        ++a;
    }
}

template <class T>
void func(const T& container)
{
    using std::begin;
    using std::end;
    func(begin(container), end(container));
}
*/

namespace StatsLibExtended
{
    template <class FWIt>
    {
    private:
        /* data */
    public:
        FWIt(FWIt a, const FWIt b);
        FWIt(const T&)
        ~FWIt();
    };
    
    FWIt::FWIt(const T& container)
    {
        using std::begin;
        using std::end;
        
    }
    FWIt::FWIt(FWIt a, const FWIt b)
    {
        //
    }
    
    FWIt::~FWIt()
    {
    }
    
}