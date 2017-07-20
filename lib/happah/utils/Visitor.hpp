// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

//pawned from http://www.artima.com/cppsource/cooperative_visitor.html

#pragma once

#include <memory>
#include <unordered_map>

template<typename Base, typename Function>
struct VTable {
     template<typename Visitable>
     void add(Function f) { m_table[typeid(Visitable).hash_code()] = f; }

protected:
     VTable() {}

     std::unordered_map<size_t, Function> m_table;
};

template<typename ReferenceVisitorVTABLE, typename Base, typename ReturnType, typename Function>
struct ReferenceVTABLE : public VTable<Base, Function> {
     ReturnType operator()(ReferenceVisitorVTABLE& visitor, Base& b) {
          auto f = this->m_table.find(typeid(b).hash_code());
          return f == this->m_table.end() ? static_cast<ReturnType>(0) : (visitor.*(f->second))(b);
     }
};

template<typename ReferenceVisitorVTABLE, typename Base, typename Function>
struct ReferenceVTABLE<ReferenceVisitorVTABLE, Base, void, Function> : public VTable<Base, Function> {
     void operator()(ReferenceVisitorVTABLE& visitor, Base& b) {
          auto f = this->m_table.find(typeid(b).hash_code());
          if(f != this->m_table.end())  (visitor.*(f->second))(b);
     }
};

template<typename SharedPointerVisitorVTABLE, typename Base, typename ReturnType, typename Function>
struct SharedPointerVTABLE : public VTable<Base, Function> {
     ReturnType operator()(SharedPointerVisitorVTABLE& visitor, std::shared_ptr<Base> b) {
          auto f = this->m_table.find(typeid(*b).hash_code());
          return f == this->m_table.end() ? static_cast<ReturnType>(0) : (visitor.*(f->second))(b);
     }
};

template<typename SharedPointerVisitorVTABLE, typename Base, typename Function>
struct SharedPointerVTABLE<SharedPointerVisitorVTABLE, Base, void, Function> : public VTable<Base, Function> {
     void operator()(SharedPointerVisitorVTABLE& visitor, std::shared_ptr<Base> b) {
          auto f = this->m_table.find(typeid(*b).hash_code());
          if(f != this->m_table.end()) (visitor.*(f->second))(b);
     }
};


template<typename Visitable, typename Base>
struct GetVisitMethodArgumentType {
     typedef Visitable TYPE;
};

template<typename Visitable, typename Base>
struct GetVisitMethodArgumentType<Visitable, const Base> {
     typedef const Visitable TYPE;
};

template<typename Visitor, typename Visitable>
struct CreateVTableBase {
     CreateVTableBase(typename Visitor::VTABLE& vtable) { vtable.template add<Visitable>(&Visitor::template thunk<Visitable, typename Visitor::INVOKER>); }
};

template<typename Visitor, typename... Visitables>
struct CreateVTable : CreateVTableBase<Visitor,Visitables>... {
     CreateVTable(typename Visitor::VTABLE& vtable) : CreateVTableBase<Visitor,Visitables>(vtable)... {}
};

template<typename VisitorImpl, typename Base, typename ReturnType, typename... Visitables>
class ReferenceVisitorVTABLE {
public:
     typedef ReturnType (ReferenceVisitorVTABLE::*FUNC)(Base&);
     typedef ReferenceVTABLE<ReferenceVisitorVTABLE, Base, ReturnType, FUNC> VTABLE;

     ReferenceVisitorVTABLE() { CreateVTable<VisitorImpl, Visitables...> createVTable(m_vtable); }

     template<typename Visitable, typename Invoker>
     ReturnType thunk(Base& b) {
          typedef typename GetVisitMethodArgumentType<Visitable, Base>::TYPE VisitableType;
          VisitorImpl& visitor = static_cast<VisitorImpl&>(*this); 
          VisitableType& visitable = static_cast<VisitableType&>(b);
          return Invoker::invoke(visitor, visitable);
     }

     ReturnType operator()(Base& b) {
          return m_vtable(*this, b);
     }

private:
     VTABLE m_vtable;
};

template<typename VisitorImpl, typename Base, typename... Visitables>
class ReferenceVisitorVTABLE<VisitorImpl, Base, void, Visitables...> {
public:
     typedef void (ReferenceVisitorVTABLE::*FUNC)(Base&);
     typedef ReferenceVTABLE<ReferenceVisitorVTABLE, Base, void, FUNC> VTABLE;

     ReferenceVisitorVTABLE() { CreateVTable<VisitorImpl, Visitables...> createVTable(m_vtable); }

     template<typename Visitable, typename Invoker>
     void thunk(Base& b) {
          typedef typename GetVisitMethodArgumentType<Visitable, Base>::TYPE VisitableType;
          VisitorImpl& visitor = static_cast<VisitorImpl&>(*this); 
          VisitableType& visitable = static_cast<VisitableType&>(b);
          Invoker::invoke(visitor, visitable);
     }

     void operator()(Base& b) {
          m_vtable(*this, b);
     }

private:
     VTABLE m_vtable;
};


template<typename VisitorImpl, typename Base, typename ReturnType, typename... Visitables>
class SharedPointerVisitorVTABLE {
public:
     typedef ReturnType (SharedPointerVisitorVTABLE::*FUNC)(std::shared_ptr<Base>);
     typedef SharedPointerVTABLE<SharedPointerVisitorVTABLE, Base, ReturnType, FUNC> VTABLE;

     SharedPointerVisitorVTABLE() { CreateVTable<VisitorImpl, Visitables...> createVTable(m_vtable); }

     template<typename Visitable, typename Invoker>
     ReturnType thunk(std::shared_ptr<Base> b) {
          typedef typename GetVisitMethodArgumentType<Visitable, Base>::TYPE VisitableType;
          VisitorImpl& visitor = static_cast<VisitorImpl&>(*this); 
          std::shared_ptr<VisitableType> visitable = std::static_pointer_cast<VisitableType>(b);
          return Invoker::invoke(visitor, visitable);
     }

     ReturnType operator()(std::shared_ptr<Base> b) {
          return m_vtable(*this, b);
     }

private:
     VTABLE m_vtable;
};

template<typename VisitorImpl, typename Base, typename... Visitables>
class SharedPointerVisitorVTABLE<VisitorImpl, Base, void, Visitables...> {
public:
     typedef void (SharedPointerVisitorVTABLE::*FUNC)(std::shared_ptr<Base>);
     typedef SharedPointerVTABLE<SharedPointerVisitorVTABLE, Base, void, FUNC> VTABLE;

     SharedPointerVisitorVTABLE() { CreateVTable<VisitorImpl, Visitables...> createVTable(m_vtable); }

     template<typename Visitable, typename Invoker>
     void thunk(std::shared_ptr<Base> b) {
          typedef typename GetVisitMethodArgumentType<Visitable, Base>::TYPE VisitableType;
          VisitorImpl& visitor = static_cast<VisitorImpl&>(*this); 
          std::shared_ptr<VisitableType> visitable = std::static_pointer_cast<VisitableType>(b);
          Invoker::invoke(visitor, visitable);
     }

     void operator()(std::shared_ptr<Base> b) {
          m_vtable(*this, b);
     }

private:
     VTABLE m_vtable;
};

#define VISIT_INVOKER(invokerName,returnType,methodName) \
struct invokerName {\
     template<typename VisitorImpl, typename Visitable>\
     static returnType invoke(VisitorImpl& visitor, Visitable& visitable) {\
          return visitor.methodName(visitable);\
     }\
};

