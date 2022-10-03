
// Copyright

#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include <queue>
#include "Point.hpp"

template <size_t N, typename ElemType>
class KDTree {
 public:
  typedef std::pair<Point<N>, ElemType> value_type;

  KDTree();

  ~KDTree();

  KDTree(const KDTree &rhs);
  KDTree &operator=(const KDTree &rhs);

  size_t dimension() const;

  size_t size() const;
  bool empty() const;

  bool contains(const Point<N> &pt) const;

  void insert(const Point<N> &pt, const ElemType &value);

  ElemType &operator[](const Point<N> &pt);

  ElemType &at(const Point<N> &pt);
  const ElemType &at(const Point<N> &pt) const;

  ElemType knn_value(const Point<N> &key, size_t k) const;

  std::vector<ElemType> knn_query(const Point<N> &key, size_t k) const;
  //void knn_queryR(const Point<N>& key, size_t k, std::priority_queue<P, std::vector<P>, std::greater<P>>& knn) const
  
 private:
  size_t dimension_;
  size_t size_;

  struct KDTNode {
      value_type v;
      KDTNode() {}
      KDTNode(const value_type& value) {
          v = value;
          children[0] = children[1] = 0;
      }
      ~KDTNode() {
          delete children[0];
          delete children[1];
      }
      KDTNode* children[2]; 
  };

  KDTNode* root;

  typename KDTNode*const* find(const Point<N>& pt) const;
  void copyRec(KDTNode*& c, const KDTNode* n);
  void knn_queryR(const Point<N>& key, size_t k, KDTNode* r, std::priority_queue<std::pair<double, ElemType>>& knn, int& min_dist, int dim) const;

};

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
    size_ = 0;
    dimension_ = N;
    root = 0;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
    delete root;
    size_ = 0;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree &rhs) {
    (*this) = rhs;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::copyRec(KDTNode*& c, const KDTNode* n) {
    if (n!=0) {
        c = new KDTNode(n->v);
        copyRec(const_cast<KDTNode*>(c->children[0]), n->children[0]);
        copyRec(const_cast<KDTNode*>(c->children[1]), n->children[1]);
    }
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>& KDTree<N, ElemType>::operator=(const KDTree& rhs) {
    copyRec(const_cast<KDTNode*>(root), rhs.root);
    size_ = rhs.size_;
    dimension_ = rhs.dimension_;
    return (*this);
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
  return dimension_;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
    return size_;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {
  return size_==0;
}

template <size_t N, typename ElemType>
typename KDTree<N, ElemType>::KDTNode*const* KDTree<N, ElemType>::find(const Point<N>& pt) const {
    int d = 0;
    auto p = &root;
    for (;
        *p && (*p)->v.first != pt;
        p = &((*p)->children[(*p)->v.first[d%N] < pt[d%N]]) , d++ //va cambiando la dimensión
        );
    //return *p != 0;
    return p;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N> &pt) const {
    auto f = find(pt);
    return *f!=0;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N> &pt, const ElemType &value) {
    if (!root) {
        std::pair< Point<N>, ElemType> aux = { pt, value };
        root = new KDTNode(aux);
        size_++;
    }
    else {
        auto f = find(pt);
        if (*f == 0) {
            std::pair< Point<N>, ElemType> aux = { pt, value };
            const_cast<KDTNode*>(*f) = new KDTNode(aux);
            size_++;
        }
        else {
            (*f)->v.first = pt;
            (*f)->v.second = value;
        }
    }
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::operator[](const Point<N> &pt) {
    auto f = find(pt);
    if (*f==0){
        ElemType e;
        insert(pt, e);
    }
    return (*f)->v.second;
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) {
    auto f = find(pt);
    if ((*f) != 0)
        return (*f)->v.second;
    throw std::out_of_range("No se encontró el valor en el árbol");
}

template <size_t N, typename ElemType>
const ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) const {
    auto f = find(pt);
    if ((*f) != 0)
        return (*f)->v.second;
    throw std::out_of_range("No se encontró el valor en el árbol");
}

template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::knn_value(const Point<N> &key, size_t k) const {
    auto v = knn_query(key, k);
    ElemType new_element;
    sort(v.begin(), v.end());
    int max = 0;
    ElemType elem = v[0];
    int cuenta = 0;
    for (int i = 1; i < v.size(); i++) {
        //std::cout << v[i] << " " << k << std::endl;
        if (v[i] == v[i - 1])
            cuenta++;
        else
            cuenta = 1;
        if (cuenta > max) {
            max = cuenta;
            elem = v[i - 1];
        }
    }
    new_element = elem;
    return new_element;
}

template <size_t N, typename ElemType>
std::vector<ElemType> KDTree<N, ElemType>::knn_query(const Point<N>& key, size_t k) const {
    typedef std::pair<double, ElemType> P;
    std::vector<ElemType> values;
    std::priority_queue<P> knn;
    int dim = 0;
    int min_dist = 1000000;
    knn_queryR(key,k,root,knn, min_dist, dim);
    while (!knn.empty()) { values.push_back(knn.top().second); knn.pop(); }
    //std::cout << "v size: " << values.size() << std::endl;
    return values;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::knn_queryR(const Point<N>& key, size_t k, KDTNode* r, std::priority_queue<std::pair<double, ElemType>> &knn, int &min_dist, int dim) const{
    typedef std::pair<double, ElemType> P;
    if (r != 0) {
        auto p = &r;
        double d = distance((*p)->v.first, key); //std::cout << d << std::endl;
        if (d < min_dist) 
            min_dist = d;
        P aux = { d,(*p)->v.second };
        knn.push(aux);
        while (knn.size() > k) knn.pop();
         
        int dir = (*p)->v.first[dim % N] < key[dim % N];
        knn_queryR(key, k, (*p)->children[dir], knn, min_dist, dim+1);

        if (fabs((*p)->v.first[dim % N] - key[dim % N]) < min_dist || knn.size() < k){ //agregar también si faltan nodos para k
          knn_queryR(key, k, (*p)->children[!dir], knn, min_dist, dim+1);
    }
   
   }

}

// TODO(me): finish the implementation of the rest of the KDTree class

#endif  // SRC_KDTREE_HPP_
