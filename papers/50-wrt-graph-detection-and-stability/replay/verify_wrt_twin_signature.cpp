#include<bits/stdc++.h>
using namespace std;
using ll=long long;
constexpr int A=4,B=4,N0=8,MAXN=11;
struct Key {
  int d,j,a,k; string monomial;
  bool operator<(Key const&o)const {return tie(d,j,a,k,monomial)<tie(o.d,o.j,o.a,o.k,o.monomial);}
};
map<Key,ll> diff;
int sign,d,a,kmax,n,m,pinmask;
int pos[MAXN], rightmost_A, forward_mask[32], forward_count;
vector<pair<int,int>> es;
vector<int> chosen;
void record(int union_mask) {
  // Only touched twin prototypes are enumerated. All other twins are inserted
  // later, in available hard-order gaps after the rightmost A block.
  for(int i=N0;i<m;i++)if(!(pinmask>>i&1)) {
    int p=pos[i];
    if (!(union_mask & ((1<<(p-1))|(1<<p)))) return;
  }
  int p=rightmost_A;
  while(p<=n && (union_mask>>p&1))p++;
  if (p==n+1)return; // A in the fixed L endpoint: no untouched twins permitted.
  int j=0;
  for(int q=p+1;q<=n;q++)if(!(union_mask>>(q-1)&1))j++;
  assert(j<=B);
  int k=__builtin_popcount((unsigned)union_mask), z=0,v[3]={};
  for(int s=union_mask;s;s&=s-1)v[z++]=s&-s;
  int perm[3]={0,1,2};string best;
  do{
    string s;
    for(int e:chosen){
      int w=0;
      for(int q=0;q<k;q++)if(v[perm[q]]&abs(e))w|=1<<q;
      s.push_back((e>0?'a':'A')+w-1);
    }
    sort(s.begin(),s.end());
    if(best.empty()||s<best)best=s;
  }while(next_permutation(perm,perm+k));
  diff[{d,j,a,k,best}]+=sign;
}
void dfs(int from,int um) {
  record(um);
  for(int i=from;i<forward_count;i++){
    int v=um|forward_mask[i];
    if(__builtin_popcount((unsigned)v)>kmax)continue;
    chosen.push_back(forward_mask[i]);dfs(i+1,v);chosen.pop_back();
  }
}
void scan_order(){
  int req=0;forward_count=0;chosen.clear();
  for(auto [v,w]:es){
    if((pinmask>>v&1) && (pinmask>>w&1))continue;
    int p=pos[v],q=pos[w];
    int mask=((1<<max(p,q))-1)^((1<<min(p,q))-1);
    if(p>q){
      req|=mask;
      if(__builtin_popcount((unsigned)req)>kmax)return;
      chosen.push_back(-mask);
    }else if(__builtin_popcount((unsigned)mask)<=kmax)forward_mask[forward_count++]=mask;
  }
  rightmost_A=*max_element(pos,pos+A);
  dfs(0,req);
}
void graph(array<int,A> rows,int s){
  sign=s;
  for(d=0;d<=3;d++){
    m=N0+d;es.clear();
    for(int i=0;i<A;i++)for(int j=0;j<B+d;j++)if(j>=B||rows[i]>>j&1)es.push_back({i,A+j});
    for(pinmask=0;pinmask<(1<<m);pinmask++){
      a=__builtin_popcount((unsigned)pinmask);if(a>3)continue;
      kmax=3-a;n=m-a;
      vector<int> freev;
      for(int i=0;i<m;i++){
        if(!(pinmask>>i&1))freev.push_back(i);
        else pos[i]=(i<A?0:n+1);
      }
      do {
        for(int i=0;i<n;i++)pos[freev[i]]=i+1;
        scan_order();
      }while(next_permutation(freev.begin(),freev.end()));
    }
  }
}
int main(){
  graph({3,7,11,12},1);
  graph({3,5,14,14},-1);
  for(auto [K,v]:diff)if(v)cout<<K.d<<' '<<K.j<<' '<<K.a<<' '<<K.k<<' '<<(K.monomial.empty()?"-":K.monomial)<<' '<<v<<std::endl;
}
