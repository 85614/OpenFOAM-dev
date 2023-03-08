#include "incompressibleFluid.H"

struct record_t
{

    using ptr_t = std::unique_ptr<void, std::function<void(void *)>>;

    template <typename _Ty>
    static ptr_t make_ptr(_Ty *ptr)
    {
        return ptr_t(ptr, [](void *p)
                     { delete static_cast<_Ty *>(p); });
    }

    std::unordered_map<std::string, ptr_t> _record;

    template <typename _Ty>
    void insert(const char *name, _Ty *_Ptr)
    {
        auto it = _record.emplace(name, make_ptr(_Ptr));
        assert(it.second);
    }

    template <typename _Ty>
    void insert_clone(const char *name, const _Ty &_Val)
    {
        auto it = _record.emplace(name, make_ptr(new _Ty(_Val)));
        assert(it.second);
    }

    template <typename _Ty>
    static auto _equal(_Ty first, _Ty second)
        -> typename std::enable_if<std::is_arithmetic<typename std::decay<_Ty>::type>::value, bool>::type
    {
        // 算术类型
        return first == second;
    }

    template <typename _Ty>
    static bool _equal(const Foam::UList<_Ty> &first, const Foam::UList<_Ty> &second)
    {
        // 列表 OpenFOAM 有实现
        return (first == second);
    }

    template <typename _Ty>
    static auto _equal(const _Ty &first, const _Ty &second)
        -> typename std::enable_if<!Foam::False<decltype(first.boundaryField)>, bool>::type
    {
        // 有 boundaryField 的
        return _equal(first.internalField()) &&
               _equal(first.boundaryField());
    }

    template <typename _Ty>
    static auto _equal(const _Ty &first, const _Ty &second)
        -> typename std::enable_if<!Foam::False<decltype(_equal(first[0], second[0]))>, bool>::type
    {
        //
        forAll(first, i)
        {
            if (!_equal(first[i], second[i]))
                return false;
        }
        return true;
    }

    template <typename _Ty>
    static bool _equal(const Foam::fvMatrix<_Ty> &first, const Foam::fvMatrix<_Ty> &second)
    {

        auto &cnfirst = const_cast<Foam::fvMatrix<_Ty> &>(first);
        auto &cnsecond = const_cast<Foam::fvMatrix<_Ty> &>(second);
        return first.hasLower() == second.hasLower() &&
               first.hasDiag() == second.hasDiag() &&
               first.hasUpper() == second.hasUpper() &&
               (!first.hasLower() || _equal(first.lower(), second.lower())) &&
               (!first.hasDiag() || _equal(first.diag(), second.diag())) &&
               (!first.hasUpper() || _equal(first.upper(), second.upper())) &&
               _equal(first.source(), second.source()) &&
               _equal(cnfirst.internalCoeffs(), cnsecond.internalCoeffs()) &&
               _equal(cnfirst.boundaryCoeffs(), cnsecond.boundaryCoeffs());
    }

    template <typename _Ty>
    bool equal(const char *name, const _Ty &_Val)
    {
        return _equal(_Val, *static_cast<_Ty *>(_record[name].get()));
    }

    template <typename _Ty>
    void assign(const char *name, const _Ty &_Val)
    {
        *static_cast<_Ty *>(_record[name].get()) = _Val;
    }

    template <typename _Ty>
    _Ty &get(const char *name)
    {
        return *static_cast<_Ty *>(_record[name].get());
    }
} record;
#define UPDATE(x) record.assign(#x, (x))
#define INSERT(x) record.insert(#x, (x))
#define EQUAL(x) record.equal(#x, (x))
#define PAIR(x) #x, (x)
