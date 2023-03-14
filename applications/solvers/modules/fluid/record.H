#include <memory>
#include <cassert>
#include <functional>
#include "solver.H"

void printCaller();

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
    bool _insert(const char *name, _Ty *_Ptr, const char *func, const char *file, int line)
    {
        Foam::Info << "init " << name << Foam::endl;
        if (func)
        {
            Foam::Info << "at " << func << " in " << file << ":" << line << Foam::endl;
        }

        auto it = _record.emplace(name, make_ptr(_Ptr));
        return it.second;
    }

    template <typename _Ty>
    void insert(const char *name, _Ty *_Ptr, const char *func, const char *file, int line)
    {
        assert(_insert(name, _Ptr, func, file, line));
    }

    template <typename _Ty>
    void insert_clone(const char *name, const _Ty &_Val, const char *func, const char *file, int line)
    {
        insert(name, new _Ty(_Val), func, file, line);
    }

    template <typename _Ty>
    void insert_clone_or_assign(const char *name, const _Ty &_Val, const char *func, const char *file, int line)
    {
        // insert(name, new _Ty(_Val));
        auto it = _record.find(name);
        if (it == _record.end())
            insert_clone(name, _Val, func, file, line);
        else
            assign(name, _Val, func, file, line);
    }

    template <typename _Ty>
    static auto _assign(_Ty &ref, const _Ty &_Val) -> typename std::enable_if<std::is_void<decltype(ref == ref)>::value>::type
    {
        ref == _Val;
        // ref = _Val 会变成 ref = std::move(_Val)!
    }

    template <typename _Ty>
    static auto _assign(_Ty &ref, const _Ty &_Val) -> typename std::enable_if<!std::is_void<decltype(ref == ref)>::value>::type
    {
        ref = _Val;
    }

    template <typename _Ty>
    void assign(const char *name, const _Ty &_Val, const char *func, const char *file, int line)
    {
        Foam::Info << "update " << name << Foam::endl;
        if (func)
        {
            Foam::Info << "at " << func << " in " << file << ":" << line << Foam::endl;
        }
        _assign(*static_cast<_Ty *>(_record[name].get()), _Val);
    }
    template <typename _Ty>
    typename std::decay<_Ty>::type &get(const char *name)
    {
        return *static_cast<typename std::decay<_Ty>::type *>(_record[name].get());
    }

    template <typename _Ty>
    _Ty &get(const char *name, const _Ty &, ...)
    {
        return get<_Ty>(name);
    }

    template <typename _Ty>
    static auto __equal(_Ty first, _Ty second)
        -> typename std::enable_if<std::is_arithmetic<typename std::decay<_Ty>::type>::value, bool>::type
    {
        // 算术类型
        return first == second;
    }

    template <typename _Ty>
    static bool __equal(const Foam::UList<_Ty> &first, const Foam::UList<_Ty> &second)
    {
        // 列表 OpenFOAM 有实现
        return (first == second);
    }

    template <class Type, template <class> class PatchField, class GeoMesh>
    static bool _equal(const Foam::GeometricField<Type, PatchField, GeoMesh> &first, const Foam::GeometricField<Type, PatchField, GeoMesh> &second)
    {
        // 有 boundaryField 的
        return _equal(first.internalField(), second.internalField()) &&
               _equal(first.boundaryField(), second.boundaryField());
    }

    template <typename _Ty>
    static auto __equal(const _Ty &first, const _Ty &second)
        -> typename std::enable_if<!Foam::False<decltype(__equal(first[0], second[0]))>, bool>::type
    {
        //
        forAll(first, i)
        {
            if (!__equal(first[i], second[i]))
                return false;
        }
        return true;
    }

    template <typename _Ty>
    static bool _equal(const Foam::tmp<_Ty> &first, const Foam::tmp<_Ty> &second)
    {
        return _equal(first(), second());
    }

    template <typename _Ty>
    static bool _equal(const _Ty &first, const _Ty &second)
    {
        return __equal(first, second);
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
    bool equal(const char *name, const _Ty &_Val, ...)
    {
        return _equal(_Val, *static_cast<_Ty *>(_record[name].get()));
    }

}; 

extern record_t record;
#define UPDATE(x) record.assign(#x, (x))
#define INSERT(x) record.insert(#x, (x))
#define EQUAL(x) record.equal(#x, (x))
#define _WHERE __FUNCTION__, __FILE__, __LINE__
#define PAIR(x) #x, (x), _WHERE
