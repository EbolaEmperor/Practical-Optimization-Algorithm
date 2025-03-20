#pragma once

#include <vector>

struct avlNode{
    int val;
    int mxval;
    int size;
    avlNode *lson, *rson;
};

class avlTree{
private:
    avlNode *root;
    void build(avlNode* &rt, const int &l, const int &r, const std::vector<int> &vals);
    void increse(avlNode* rt, const int &k, const int &x);
    void setzero(avlNode* rt, const int &k);
    int getmaxID(avlNode* rt) const;

public:
    avlTree();
    void build(const std::vector<int> &vals);
    void increse(const int &k, const int &x);
    void setzero(const int &k);
    int getmaxID() const;
};


#include <vector>
using std::vector;

avlTree::avlTree(){
    root = nullptr;
}

void upmax(int &x, const int &y){
    if(y > x) x = y;
}

void avlTree::build(avlNode* &rt, const int &l, const int &r, const vector<int> &vals){
    if(rt == nullptr){
        rt = new avlNode;
        rt->lson = rt->rson = nullptr;
    }
    int mid = (l+r)/2;
    rt->mxval = rt->val = vals[mid];
    if(l < mid){
        build(rt->lson, l, mid-1, vals);
        upmax(rt->mxval, rt->lson->mxval);
    }
    if(r > mid){
        build(rt->rson, mid+1, r, vals);
        upmax(rt->mxval, rt->rson->mxval);
    }
    rt->size = r-l+1;
}

void avlTree::increse(avlNode* rt, const int &k, const int &x){
    if(!rt->lson && k==0 || rt->lson && k == rt->lson->size){
        rt->val += x;
        upmax(rt->mxval, rt->val);
        return;
    }
    if(!rt->lson || k > rt->lson->size){
        increse(rt->rson, rt->lson ? k-(rt->lson->size+1) : k-1, x);
        upmax(rt->mxval, rt->rson->mxval);
    } else {
        increse(rt->lson, k, x);
        upmax(rt->mxval, rt->lson->mxval);
    }
}

void avlTree::setzero(avlNode* rt, const int &k){
    if(!rt->lson && k==0 || rt->lson && k == rt->lson->size){
        rt->val = -1;
    } else if(!rt->lson || k > rt->lson->size){
        setzero(rt->rson, rt->lson ? k-(rt->lson->size+1) : k-1);
    } else {
        setzero(rt->lson, k);
    }
    rt->mxval = rt->val;
    if(rt->lson) upmax(rt->mxval, rt->lson->mxval);
    if(rt->rson) upmax(rt->mxval, rt->rson->mxval);
}

int avlTree::getmaxID(avlNode* rt) const{
    if(rt->lson && rt->lson->mxval == rt->mxval){
        return getmaxID(rt->lson);
    } else if(rt->val == rt->mxval){
        return rt->lson ? rt->lson->size : 0;
    } else {
        return getmaxID(rt->rson) + (rt->lson ? rt->lson->size : 0) + 1;
    }
}

void avlTree::build(const vector<int> &vals){
    build(root, 0, vals.size()-1, vals);
}

void avlTree::increse(const int &k, const int &x){
    increse(root, k, x);
}

void avlTree::setzero(const int &k){
    setzero(root, k);
}

int avlTree::getmaxID() const{
    return getmaxID(root);
}
