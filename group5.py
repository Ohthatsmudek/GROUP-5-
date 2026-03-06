# group5.py - Run: python group5.py
import math
try:
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    HAS_PLOT = True
except ImportError:
    HAS_PLOT = False

def mat_inv2x2(M):
    a,b=M[0][0],M[0][1]; c,d=M[1][0],M[1][1]
    det=a*d-b*c
    if abs(det)<1e-15: raise ValueError("Singular matrix.")
    return [[d/det,-b/det],[-c/det,a/det]]

def least_squares_2x2(A,b_vec):
    n=len(A); AtA=[[0.0,0.0],[0.0,0.0]]; Atb=[0.0,0.0]
    for i in range(n):
        for r in range(2):
            for c in range(2): AtA[r][c]+=A[i][r]*A[i][c]
            Atb[r]+=A[i][r]*b_vec[i]
    inv=mat_inv2x2(AtA)
    return inv[0][0]*Atb[0]+inv[0][1]*Atb[1], inv[1][0]*Atb[0]+inv[1][1]*Atb[1]

def get_float(p):
    while True:
        try: return float(input(p))
        except ValueError: print("  >> Enter a valid number.")

def get_int(p,minimum=1):
    while True:
        try:
            v=int(input(p))
            if v<minimum: print(f"  >> Minimum is {minimum}.")
            else: return v
        except ValueError: print("  >> Enter a whole number.")

def get_dms(p):
    while True:
        try:
            parts=input(f"{p} (D M S): ").strip().split()
            if len(parts)!=3: raise ValueError
            d,m,s=map(float,parts)
            return d+m/60.0+s/3600.0
        except ValueError: print("  >> Enter three numbers e.g.: 47 15 30")

def dms_str(deg):
    deg=abs(deg); d=int(deg); m=int((deg-d)*60)
    s=round(((deg-d)*60-m)*60,2)
    return f"{d}\u00b0 {m}' {s}\""

def sep(title=""):
    line="="*72
    if title: print(f"\n{line}\n  {title}\n{line}")
    else: print(line)

def tbl(headers,rows,widths=None):
    if widths is None:
        widths=[len(h) for h in headers]
        for row in rows:
            for i,cell in enumerate(row): widths[i]=max(widths[i],len(str(cell)))
    fmt="  ".join(f"{{:<{w}}}" for w in widths)
    div="  ".join("-"*w for w in widths)
    print(fmt.format(*[str(h) for h in headers])); print(div)
    for row in rows: print(fmt.format(*[str(c) for c in row]))
    print(div)

def compute_provisional(name_P,stations,idx_A,idx_B):
    s1=stations[idx_A]; s2=stations[idx_B]
    if s1[3]>=s2[3]: name_A,N_A,E_A,brg_A=s1; name_B,N_B,E_B,brg_B=s2
    else: name_A,N_A,E_A,brg_A=s2; name_B,N_B,E_B,brg_B=s1
    eps=1e-9; brg_A_rad=math.radians(brg_A); brg_B_rad=math.radians(brg_B)
    if abs(math.cos(brg_A_rad))<eps or abs(math.cos(brg_B_rad))<eps:
        print("\nERROR: One bearing near 90/270. Choose different stations."); return None,None
    t_A=math.tan(brg_A_rad); t_B=math.tan(brg_B_rad); denom=t_B-t_A
    if abs(denom)<eps: print("\nERROR: Bearings parallel."); return None,None
    rhs_A=E_A-N_A*t_A; rhs_B=E_B-N_B*t_B
    N_P=(rhs_A-rhs_B)/denom; E_P=E_A+(N_P-N_A)*t_A; E_P_check=E_B+(N_P-N_B)*t_B
    dist_AP=math.hypot(N_P-N_A,E_P-E_A); dist_BP=math.hypot(N_P-N_B,E_P-E_B)
    dN_AP=N_P-N_A; dE_AP=E_P-E_A; dN_BP=N_P-N_B; dE_BP=E_P-E_B
    dN_AB=N_B-N_A; dE_AB=E_B-E_A; delta_alpha=brg_A-brg_B
    sep("COMPUTATION TABLE  -  BEARING METHOD")
    W=[30,30,30,28]
    rows=[
        [f"{name_P}  (unknown station)",f"Nc = {N_A:.2f} + ({dN_AP:+.2f})",f"Ec = {E_A:.2f} + ({dE_AP:+.2f})",f"Dist {name_A}-{name_P} = {dist_AP:.2f} m"],
        ["",f"  dN = {dist_AP:.2f} x cos({brg_A:.2f})",f"  dE = {dist_AP:.2f} x sin({brg_A:.2f})",""],
        ["",f"     = {dN_AP:+.2f}",f"     = {dE_AP:+.2f}",""],["","","",""],
        [f"{name_A}  (larger bearing)",f"N = {N_A:.2f}",f"E = {E_A:.2f}",f"a_{name_A}{name_P} = {dms_str(brg_A)}"],
        ["",f"  dN({name_A}{name_B}) = {dN_AB:+.2f}",f"  dE({name_A}{name_B}) = {dE_AB:+.2f}",f"  da  =  {dms_str(brg_A)}"],
        ["","","",f"        -  {dms_str(brg_B)}"],["","","",f"        =  {dms_str(delta_alpha)}"],["","","",""],
        [f"{name_B}  (lesser bearing)",f"N = {N_B:.2f}",f"E = {E_B:.2f}",f"a_{name_B}{name_P} = {dms_str(brg_B)}"],
        ["",f"  dN = {dist_BP:.2f} x cos({brg_B:.2f})",f"  dE = {dist_BP:.2f} x sin({brg_B:.2f})",f"Dist {name_B}-{name_P} = {dist_BP:.2f} m"],
        ["",f"     = {dN_BP:+.2f}",f"     = {dE_BP:+.2f}",""],["","","",""],
        [f"{name_P}  (check)",f"Nc = {N_B:.2f} + ({dN_BP:+.2f})",f"Ec = {E_B:.2f} + ({dE_BP:+.2f})",""],
        ["",f"   = {N_P:.2f}",f"   = {E_P_check:.2f}",""],
    ]
    tbl(["Station","Northings","Eastings","Distance / Bearing"],rows,widths=W)
    ok=abs(E_P-E_P_check)<0.002
    print(f"\n  Check: N = {N_P:.2f},  E = {E_P_check:.2f}  {'[CHECK OK]' if ok else '[WARNING]'}")
    print(f"\n  Provisional {name_P}:  N = {N_P:.2f},  E = {E_P:.2f}")
    return N_P,E_P

def print_cut_computation(name_stn,N_stn,E_stn,name_P,N_P,E_P,brg_deg):
    brg_rad=math.radians(brg_deg); sin_b=math.sin(brg_rad); cos_b=math.cos(brg_rad); eps=1e-9
    dN=N_P-N_stn; dE=E_P-E_stn
    cot_a=cos_b/sin_b if abs(sin_b)>eps else float('inf')
    tan_a=sin_b/cos_b if abs(cos_b)>eps else float('inf')
    dN_prime=cot_a*dE; dE_prime=tan_a*dN
    cut_N=dN_prime-dN; cut_E=dE_prime-dE
    S1=dE/sin_b if abs(sin_b)>eps else float('inf')
    S2=dN/cos_b if abs(cos_b)>eps else float('inf')
    line="-"*72
    print(f"\n{line}"); print(f"  Cut Computation for {name_stn}"); print(line)
    print(f"  {name_stn:<6}  {N_stn:.2f}      {E_stn:.2f}")
    print(f"  {name_P:<6}  {N_P:.2f}      {E_P:.2f}"); print()
    print(f"  \u0394N{name_stn.lower()}{name_P.lower()} = {dN:.2f}                  \u0394E{name_stn.lower()}{name_P.lower()} = {dE:.2f}"); print()
    print(f"  Cut N")
    print(f"  \u0394N'{name_stn.lower()}{name_P.lower()} = {dN_prime:.2f}            cot\u03b1{name_stn.lower()}{name_P.lower()} = {cot_a}")
    print(f"  cut = \u0394N'{name_stn.lower()}{name_P.lower()} - \u0394N{name_stn.lower()}{name_P.lower()} = {dN_prime:.2f} - ({dN:.2f})")
    print(f"  cut = {cut_N:.2f}")
    print(f"  S\u2081 = \u0394E cosec\u03b1{name_stn.lower()}{name_P.lower()} = {S1:.2f} m"); print()
    print(f"  Cut E")
    print(f"  \u0394E'{name_stn.lower()}{name_P.lower()} = {dE_prime:.2f}            tan\u03b1{name_stn.lower()}{name_P.lower()} = {tan_a}")
    print(f"  cut = \u0394E'{name_stn.lower()}{name_P.lower()} - \u0394E{name_stn.lower()}{name_P.lower()} = {dE_prime:.2f} - ({dE:.2f})")
    print(f"  cut = {cut_E:.2f}")
    print(f"  S\u2082 = \u0394N sec\u03b1{name_stn.lower()}{name_P.lower()} = {S2:.2f} m"); print(line)
    return cut_N,cut_E,S1,S2

def least_squares_adjust(cut_stations,N_P,E_P):
    A,b_vec=[],[]
    for _,N_C,E_C,brg in cut_stations:
        sin_b=math.sin(math.radians(brg)); cos_b=math.cos(math.radians(brg))
        A.append([-sin_b,cos_b]); b_vec.append(sin_b*(N_P-N_C)-cos_b*(E_P-E_C))
    return least_squares_2x2(A,b_vec)

def plot_cut_graph(stn_names,cut_N_vals,cut_E_vals,s_vals,name_P,N_P,E_P):
    if not HAS_PLOT:
        print("\n  [Chart skipped - matplotlib not installed]"); return None,None
    max_s=max(abs(s) for s in s_vals)
    scales=[abs(s)/max_s for s in s_vals]
    print(f"\n  Scale factors (S / max S):")
    for name,s,sc in zip(stn_names,s_vals,scales):
        print(f"    Station {name}: S = {s:.2f}  |  scale = {sc:.6f}")
    positions_N=list(cut_N_vals); positions_E=list(cut_E_vals)
    step=0.0001; max_iter=100000; tolerance=1e-6
    for _ in range(max_iter):
        total_scale=sum(scales)
        mean_N=sum(positions_N[i]*scales[i] for i in range(len(scales)))/total_scale
        mean_E=sum(positions_E[i]*scales[i] for i in range(len(scales)))/total_scale
        new_N=[positions_N[i]+(mean_N-positions_N[i])*step*scales[i] for i in range(len(scales))]
        new_E=[positions_E[i]+(mean_E-positions_E[i])*step*scales[i] for i in range(len(scales))]
        max_spread=max(max(abs(new_N[i]-mean_N) for i in range(len(scales))),
                       max(abs(new_E[i]-mean_E) for i in range(len(scales))))
        positions_N=new_N; positions_E=new_E
        if max_spread<tolerance: break
    total_scale=sum(scales)
    conv_N=sum(positions_N[i]*scales[i] for i in range(len(scales)))/total_scale
    conv_E=sum(positions_E[i]*scales[i] for i in range(len(scales)))/total_scale
    print(f"\n  Convergence point:")
    print(f"    Cut N at convergence = {conv_N:.2f}")
    print(f"    Cut E at convergence = {conv_E:.2f}")
    all_N=list(cut_N_vals)+[conv_N]; all_E=list(cut_E_vals)+[conv_E]
    range_N=max(all_N)-min(all_N); range_E=max(all_E)-min(all_E)
    pad_N=max(range_N*0.4,0.5); pad_E=max(range_E*0.4,0.5)
    y_min=min(all_N)-pad_N; y_max=max(all_N)+pad_N
    x_min=min(all_E)-pad_E; x_max=max(all_E)+pad_E
    def nice_step(rng):
        raw=rng/8
        mag=10**math.floor(math.log10(raw)) if raw>0 else 0.1
        for f in [1,2,2.5,5,10]:
            if raw<=f*mag: return f*mag
        return mag*10
    x_step=nice_step(x_max-x_min); y_step=nice_step(y_max-y_min)
    # GROUP FIVE: Light cream theme — purple, teal, pink
    colors=['#FF6B35','#FF0000','#FFD700','#FF8C00','#FF4500','#FFA500']
    fig,ax=plt.subplots(figsize=(11,9))
    ax.set_facecolor('#0a0f2e')
    fig.patch.set_facecolor('#050a1a')
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(x_step/5))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(y_step/5))
    ax.grid(which='minor',color='#1a3a6e',linewidth=0.4,linestyle='-')
    ax.xaxis.set_major_locator(ticker.MultipleLocator(x_step))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(y_step))
    ax.grid(which='major',color='#2a5caa',linewidth=0.9,linestyle='-')
    ax.axhline(y=0,color='#7ecfff',linewidth=1.2,zorder=3)
    ax.axvline(x=0,color='#7ecfff',linewidth=1.2,zorder=3)
    ax.tick_params(colors='white'); ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white'); ax.title.set_color('white')
    for spine in ax.spines.values(): spine.set_edgecolor('#2a5caa')
    for i,(sn,cN,cE) in enumerate(zip(stn_names,cut_N_vals,cut_E_vals)):
        color=colors[i%len(colors)]
        ax.plot([cE,conv_E],[cN,conv_N],color=color,linewidth=1.8,linestyle='-',zorder=4,label=f"Station {sn}")
        ax.plot(cE,cN,'o',color=color,markersize=6,zorder=6)
        ax.annotate(f"Station {sn}  ({cE:.2f}, {cN:.2f})",(cE,cN),
                    textcoords="offset points",xytext=(7,5),fontsize=8,color=color,fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.25',facecolor='#0a0f2e',alpha=0.9,edgecolor=color))
    ax.plot(conv_E,conv_N,'x',color='white',markersize=8,markeredgewidth=2.5,zorder=10,label='Convergence (P final)')
    ax.axhline(y=conv_N,color='#aaaaaa',linewidth=0.7,linestyle=':',zorder=3)
    ax.axvline(x=conv_E,color='#aaaaaa',linewidth=0.7,linestyle=':',zorder=3)
    ax.annotate(f"  P final\n  Cut N = {conv_N:.2f}\n  Cut E = {conv_E:.2f}",
                (conv_E,conv_N),textcoords="offset points",xytext=(10,-28),fontsize=9,fontweight='bold',color='white',
                bbox=dict(boxstyle='round,pad=0.4',facecolor='#1a3a6e',alpha=0.95,edgecolor='white'))
    ax.set_xlim(x_min,x_max); ax.set_ylim(y_min,y_max)
    ax.set_xlabel("Cut E  (m)",fontsize=11,fontweight='bold')
    ax.set_ylabel("Cut N  (m)",fontsize=11,fontweight='bold')
    ax.set_title(f"Cut N vs Cut E  -  Unknown Station {name_P}",fontsize=13,fontweight='bold',pad=12)
    ax.legend(fontsize=9,loc='best',framealpha=0.7,facecolor='#0a0f2e',labelcolor='white')
    ax.set_aspect('equal','box')
    plt.tight_layout(); plt.show()
    return conv_N,conv_E

def finish():
    print()
    surveyor=input("  Name of Surveyor      : ").strip()
    date=input("  Date of Computation  : ").strip()
    print(f"\n  Name of Surveyor      : {surveyor}")
    print(f"  Date of Computation  : {date}")
    print("\n"+"="*72)
    msg="END OF SEMI GRAPHICAL FIX FOR COMPUTATION OF FINAL COORDINATES"
    print(msg.center(72))
    print("="*72)

def main():
    sep()
    print("  WELCOME - SEMI GRAPHICAL METHOD DEVELOPED BY GROUP FIVE".center(72))
    print("  GROUP FIVE SURVEY SYSTEM".center(72))
    sep()
    print()
    name_P=input("  Enter the name of the unknown point to be fixed: ").strip() or "P"
    sep("KNOWN STATIONS")
    print("  We will compute provisional coordinates of the unknown station")
    print("  from intersection (bearing method).\n")
    n_total=get_int("  How many control stations do you have? (min 2): ",minimum=2)
    stations=[]
    for i in range(n_total):
        print(f"\n  Station {i+1}:")
        stn_name=input(f"    Control station name: ").strip() or f"Stn{i+1}"
        N_C=get_float(f"    Enter Northing coordinate of {stn_name}: ")
        E_C=get_float(f"    Enter Easting coordinate of {stn_name}: ")
        brg=get_dms(f"    Enter observed bearing from {stn_name} to unknown point {name_P}")
        stations.append((stn_name,N_C,E_C,brg))
    sep("SELECT STATIONS FOR PROVISIONAL COORDINATES")
    print(f"\n  Available stations:\n")
    for i,(sn,N,E,b) in enumerate(stations):
        print(f"    [{i+1}]  {sn:<14}  N = {N:.2f}    E = {E:.2f}    Bearing = {dms_str(b)}")
    print(f"\n  Choose two stations to compute provisional coordinates of {name_P}.")
    print(f"  The remaining stations will be used for cut adjustment.\n")
    while True:
        idx_A=get_int("  Select first station for intersection:  ",minimum=1)-1
        idx_B=get_int("  Select second station for intersection: ",minimum=1)-1
        if idx_A==idx_B: print("  >> Pick two DIFFERENT stations.")
        elif idx_A>=n_total or idx_B>=n_total: print(f"  >> Numbers must be 1 to {n_total}.")
        else: break
    sep("COMPUTATION OF PROVISIONAL COORDINATES (Bearing Method)")
    print(f"\n  Unknown station : {name_P}")
    print(f"  Using stations  : {stations[idx_A][0]}  and  {stations[idx_B][0]}\n")
    N_P,E_P=compute_provisional(name_P,stations,idx_A,idx_B)
    if N_P is None: finish(); return
    cut_stations=[s for i,s in enumerate(stations) if i not in (idx_A,idx_B)]
    if len(cut_stations)<2:
        print(f"\n  NOTE: Only {len(cut_stations)} station(s) left (need at least 2).")
        print(f"  Final coordinates of {name_P}:  N = {N_P:.2f},  E = {E_P:.2f}")
        finish(); return
    sep("COORDINATE CUTS")
    cuts=[]; cut_N_vals=[]; cut_E_vals=[]; s1_vals=[]; s2_vals=[]; stn_names=[]
    for sn,N_C,E_C,brg in cut_stations:
        cut_N,cut_E,S1,S2=print_cut_computation(sn,N_C,E_C,name_P,N_P,E_P,brg)
        cuts.append((cut_N,cut_E,S1)); cut_N_vals.append(cut_N); cut_E_vals.append(cut_E)
        s1_vals.append(S1); s2_vals.append(S2); stn_names.append(sn)
    try:
        cN,cE=least_squares_adjust(cut_stations,N_P,E_P)
    except ValueError as e:
        print(f"\n  Adjustment error: {e}"); finish(); return
    adj_N=N_P+cN; adj_E=E_P+cE
    sep("ADJUSTMENT RESULTS")
    print(f"  Provisional : N = {N_P:.2f},   E = {E_P:.2f}")
    print(f"  Corrections : dN = {cN:+.2f},   dE = {cE:+.2f}")
    print(f"  Adjusted    : N = {adj_N:.2f},   E = {adj_E:.2f}")
    sep("RESIDUALS")
    res_rows=[]
    for i,(sn,N_C,E_C,brg) in enumerate(cut_stations):
        sin_b=math.sin(math.radians(brg)); cos_b=math.cos(math.radians(brg))
        res=sin_b*(adj_N-N_C)-cos_b*(adj_E-E_C)
        res_rows.append([f"Stn {i+1}  ({sn})",dms_str(brg),f"{res:+.2f} m"])
    tbl(["Station","Bearing","Residual"],res_rows)
    sep("DISTANCES TO ADJUSTED POINT")
    for sn,N_C,E_C,_ in cut_stations:
        dist=math.hypot(adj_N-N_C,adj_E-E_C)
        print(f"  {sn:<14}  N = {N_C:.2f}   E = {E_C:.2f}   Dist to {name_P} = {dist:.2f} m")
    if HAS_PLOT:
        sep("CUT N vs CUT E GRAPH")
        print("\n  S1 values per station:")
        for sn,s in zip(stn_names,s1_vals): print(f"    Station {sn}: S1 = {s:.2f} m")
        print("\n  S2 values per station:")
        for sn,s in zip(stn_names,s2_vals): print(f"    Station {sn}: S2 = {s:.2f} m")
        while True:
            choice=input("\n  Which distance will you use for scaling? (S1 or S2 - enter 1 or 2): ").strip()
            if choice=='1': s_vals_chosen=s1_vals; print("  Using S1."); break
            elif choice=='2': s_vals_chosen=s2_vals; print("  Using S2."); break
            else: print("  >> Enter 1 or 2.")
        conv_N,conv_E=plot_cut_graph(stn_names,cut_N_vals,cut_E_vals,s_vals_chosen,name_P,N_P,E_P)
        if conv_N is not None:
            P_final_N=N_P+conv_N; P_final_E=E_P+conv_E
            sep("FINAL COORDINATES")
            print(f"  Provisional {name_P} : N = {N_P:.2f},   E = {E_P:.2f}")
            print(f"  Cut N (graph)        : {conv_N:.2f}")
            print(f"  Cut E (graph)        : {conv_E:.2f}")
            print(f"\n  P final {name_P}  =  ( N = {P_final_N:.2f},   E = {P_final_E:.2f} )")
    else:
        print(f"\n  [DONE]  {name_P}  =  ( N = {adj_N:.2f},  E = {adj_E:.2f} )")
    finish()

if __name__ == "__main__":
    main()
