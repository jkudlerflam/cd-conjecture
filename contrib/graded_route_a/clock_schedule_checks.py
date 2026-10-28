"""Enumerate the C8 binary-clock path and verify its edge-use budget."""
import json
from pathlib import Path


def check(n):
    word=['I']+['Z']*(n-1)
    seen=set()
    carry=[0]*(n-1)
    returning=[0]*(n-1)
    init=[0]*(n-1)
    turns=[0]*n
    while True:
        state=tuple(word)
        assert state not in seen
        seen.add(state)
        heads=[i for i,s in enumerate(word) if s in ('I','C0','C1','R0','R1')]
        assert len(heads)==1
        i=heads[0]
        h=word[i]
        if h=='I':
            if i==n-1:
                word[i]='R0'
            else:
                assert word[i+1]=='Z'
                word[i],word[i+1]='L','I'
                init[i]+=1
        elif h=='C0':
            word[i]='R1'
            turns[i]+=1
        elif h=='C1':
            if i==n-1:
                break
            digit=word[i+1]
            assert digit in ('0','1')
            word[i],word[i+1]='L','C'+digit
            carry[i]+=1
        elif i==0:
            word[i]='C'+h[1]
        else:
            assert word[i-1]=='L'
            word[i-1],word[i]='R0',h[1]
            returning[i-1]+=1
    expected=[2**(n-j) for j in range(1,n)]
    assert len(seen)==n+4*(2**n-1)
    assert carry==returning==expected
    assert init==[1]*(n-1)
    assert turns==[2**(n-i) for i in range(1,n+1)]
    return dict(n=n,vertices=len(seen),carry_by_edge=carry,
                return_by_edge=returning,turns_by_site=turns,
                final_word=word)


if __name__=='__main__':
    out={'scope':'Gate scheduling only; no universality claim.',
         'cases':[check(n) for n in range(1,11)]}
    Path(__file__).with_name('clock_schedule_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
