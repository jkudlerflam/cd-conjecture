"""Independent physical nine-letter Hermitian nearest-neighbour matrices."""
import itertools
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

LETTERS=['L0','L1','R0','R1','DA','DB','E','MA','MB']
CHARGE=np.array([0,2,0,2,-1,1,0,-1,1])
MARKER=np.array([0,0,0,0,0,0,1,1,1])


def allowed(a,b):
    if a=='L0' or a in ('E','MA','MB'):return b!='R1'
    if a=='L1':return b in ('L0','L1','R0','DA','DB')
    return b in ('R0','R1','E')


def parts():
    ix={x:i for i,x in enumerate(LETTERS)};d=9
    matrices={name:np.zeros((d*d,)*2) for name in ('walk','convert','port','swap')}
    def add(name,out,inn,value):
        row=ix[out[0]]*d+ix[out[1]];col=ix[inn[0]]*d+ix[inn[1]]
        matrices[name][row,col]+=value;matrices[name][col,row]+=value
    for a in (0,1):
        add('walk',(f'L{a}','R0'),('R0',f'R{a}'),-1)
        for species in ('A','B'):
            add('walk',(f'L{a}',f'D{species}'),(f'D{species}',f'R{a}'),-((-1)**a if species=='A' else 1))
    add('walk',('L0','E'),('R0','E'),-1)
    for species in ('A','B'):add('walk',('L0',f'M{species}'),(f'D{species}','E'),-1)
    add('convert',('DB','R0'),('DA','R1'),1)
    add('port',('MA','DB'),('E','R0'),-1)
    add('port',('MB','DA'),('E','R0'),1)
    for side in ('L','R'):add('swap',(side+'1',side+'0'),(side+'0',side+'1'),1)
    matrices['penalty']=np.diag([int(not allowed(a,b)) for a in LETTERS for b in LETTERS])
    return matrices


def cell_words(n):
    s=n+2;words=[]
    for species in ('V','A','B'):
        for j in range(s):
            for bits in range(1<<n):
                a=[(bits>>i)&1 for i in range(n)]
                if j==s-1:
                    word=[f'L{x}' for x in a]+['L0','E' if species=='V' else 'M'+species]
                else:
                    token='R0' if species=='V' else 'D'+species
                    word=[f'L{x}' for x in a[:j]]+[token]+[f'R{x}' for x in a[j:]]+['E']
                words.append(tuple(word))
    assert len(set(words))==3*s*(1<<n)
    return words


def cell(n):
    words=cell_words(n);index={x:i for i,x in enumerate(words)};ix={x:i for i,x in enumerate(LETTERS)}
    out={};d=9
    for name,local in parts().items():
        if name=='port':continue
        row=[];col=[];value=[]
        for k,w in enumerate(words):
            for j in range(len(w)-1):
                source=ix[w[j]]*d+ix[w[j+1]]
                for dest in np.flatnonzero(local[:,source]):
                    z=w[:j]+(LETTERS[dest//d],LETTERS[dest%d])+w[j+2:]
                    if z in index:row.append(index[z]);col.append(k);value.append(local[dest,source])
        out[name]=sparse.csr_matrix((value,(row,col)),shape=(len(words),)*2)
    for name,func in [('number',lambda x:x in ('L1','R1')),('left',lambda x:x=='L1')]:
        out[name]=sparse.diags([sum(func(x) for x in w) for w in words],dtype=float,format='csr')
    s=n+2;d=1<<n;v=np.sqrt(2/(s+1))*np.sin(np.pi*np.arange(1,s+1)/(s+1))
    G=np.zeros((len(words),d))
    for j in range(s):G[j*d:(j+1)*d]=v[j]*np.eye(d)
    out.update(words=words,G=G,Ev=-2*np.cos(np.pi/(s+1)),alpha=v[0])
    return out


def pair_source(n,cell_data):
    """Literal port applied to the exact vacuum codes, in charged-pair basis."""
    s=n+2;d=1<<n;words=cell_data['words'];charged=words[s*d:];index={w:i for i,w in enumerate(charged)}
    ix={x:i for i,x in enumerate(LETTERS)};V=parts()['port'];Nc=len(charged)
    S=np.zeros((Nc*Nc,d*d))
    # Vacuum source requires the left head at its terminal and right at its start.
    for a in range(d):
        for b in range(d):
            left=words[(s-1)*d+a];right=words[b]
            source=ix[left[-1]]*9+ix[right[0]]
            for dest in np.flatnonzero(V[:,source]):
                l=left[:-1]+(LETTERS[dest//9],);r=(LETTERS[dest%9],)+right[1:]
                # Logical Fock index uses the left cell as the low bits.
                S[index[l]*Nc+index[r],a+d*b]+=V[dest,source]
    return S


def literal_coefficients(n,mass=2,max_order=4):
    c=cell(n);s=n+2;d=1<<n;start=s*d
    charged_walk=c['walk'][start:,start:];charged_C=c['convert'][start:,start:]
    Nc=charged_walk.shape[0]
    A=sparse.kronsum(charged_walk+mass*sparse.eye(Nc),charged_walk+mass*sparse.eye(Nc),format='csc')
    C=sparse.kronsum(charged_C,charged_C,format='csc')
    S=pair_source(n,c);lu=splu(A);X=lu.solve(S);out={0:-S.T@X}
    for order in range(1,max_order+1):
        X=lu.solve(C@X)
        if order%2==0:out[order]=-S.T@X
        else:assert np.max(abs(S.T@X))<1e-12
    return out,dict(cell_dimension=len(c['words']),pair_excited_dimension=A.shape[0],code_dimension=d*d)
