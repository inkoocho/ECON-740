# Econ 777 

Matlab and Julia code repositories for ECON 777 class

## 1. Gresham folder 
The Matlab and Julia codes are based on [Cho and Kasa (2017, AER)](https://www.aeaweb.org/articles?id=10.1257/aer.20160665).



## 2. Monthly stock data (crsp_m.csv)
```
/* Selected variables from the CRSP monthly data file (crsp.msf)  */
%let msfvars = prc ret shrout cfacpr cfacshr;

/* Selected variables from the CRSP monthly event file (crsp.mse) */
%let msevars = ncusip exchcd shrcd ;

/* This procedure creates a Monthly CRSP dataset named “CRSP_M”   */
%crspmerge(s=m,start=&begdate,end=&enddate,sfvars=&msfvars,sevars=&msevars,filters=&sfilter);

data crsp_m; format qdate YYMMDD10.; format mdate YYMMDD10.;
set crsp_m;
	WHERE shrout >0;
	qdate = INTNX('QTR',date,0,'E');
	mdate = INTNX("MONTH",date,0,"E");
	price = abs(prc)/cfacpr;
	shares_out = shrout*cfacshr*1000;
	if shares_out<=0 then shares_out=.;
	market_cap = price*shares_out/1000000;
	label price = "Price at Period End, Adjusted";
	label shares_out = "Total Shares Outstanding, Adjusted";
	label market_cap = "Market Capitalization, x$1m";
	drop prc cfacpr shrout shrcd;
run;
```

## 3. Quarterly 13F holdings data (13f_holding.csv & 13f_asset.csv)
```
libname tfn1 "/wrds/tfn/sasdata/s34";
libname tfn2 "/wrds/tfn/sasdata/s12";
libname mfl "/wrds/mfl/sasdata";
libname ff  "/wrds/ff/sasdata";

/* Merge TR-13f S34type1 and S34type3 Sets */
/* First, Keep First Vintage with Holdings Data for Each RDATE-MGRNO Combinations */
proc sql;
  create table First_Vint
  as select distinct rdate, fdate, mgrno, mgrname
  from tfn1.s34type1
  group by mgrno, rdate
  having fdate=min(fdate)
  order by mgrno, rdate;
quit;
 
/* Marker for First and Last Quarters of Reporting & Reporting Gaps */
/* Excercise Helpful Mostly For Clean Time-Series Analysis */
/* Idea: Compute Trades only for Institutions who File Both Q(t) and Q(t-1) 13F Reports */       
data First_Vint;
  set First_Vint;
  by mgrno rdate;
  length First_Report 3;
  First_Report = (first.mgrno or intck("QTR",lag(rdate),rdate)>1);
run;
 
/* Last Report by Institutional Manager, or Missing 13F Reports in the Next Quarter(s) */
proc sort data=First_Vint nodupkey; by mgrno descending rdate; run;
data First_Vint;
  set First_Vint;
  by mgrno descending rdate;
  length Last_Report 3;
  Last_Report = (first.mgrno or intck("QTR",rdate,lag(rdate))>1);
  if ("&begdate"d <= rdate <="&enddate"d);
run;
proc sort data=First_Vint; by fdate mgrno; run;
 
/* Extract Holdings and Adjust Shares Held */
/* FDATE -Vintage Date- is used in Shares' Adjustment */
data Holdings_v1 / view=Holdings_v1;
  merge First_Vint(in=a)
    tfn1.s34type3(in=b drop=sole shared no);
  by fdate mgrno;
  if a and b and shares>0;
run;
 
/* Map TR-13F's Historical CUSIP to CRSP Unique Identifier PERMNO */
/* Keep Securities in CRSP Common Stock Universe */
proc sql;
  create view Holdings_v2 as
  select distinct a.rdate, a.fdate, a.mgrno, a.mgrname, a.type,
          a.first_report, a.last_report, b.permno, a.shares
  from Holdings_v1 as a,
    (select distinct ncusip, permno from crsp.msenames
      where not missing(ncusip)) as b
      where a.cusip=b.ncusip;
quit;
  
/* Adjust Shares using CRSP Adjustment Factors aligned at Vintage Dates */
proc sql;
  create table Holding as
  select a.rdate, a.mgrno, a.mgrname, a.type, a.first_report, a.last_report, a.permno,
    INTNX('QTR', rdate, 0, 'E') AS date format YYMMDD10.,
    a.shares*b.cfacshr as shares_adj label = "Adjusted Shares Held"
  from Holdings_v2 as a, crsp_m as b
  where a.permno=b.permno and a.fdate = b.qdate;
quit;
 
/* Sanity Checks for Duplicates - Ultimately, Should be 0 Duplicates */
/* If No Errors, then Duplicates can be due to 2 historical CUSIPs   */
/* (Separate Holdings by Same Manager) mapping to the same permno */
/* Total # of cases should around 100 (->Get Sum or Drop Dups) */
proc sort data=Holding nodupkey thread; by mgrno permno rdate; where shares_adj>0; run;

/* Calculate Total Assets per Institution Every Quarter */
proc sql;
create table Asset
as select distinct a.mgrno, a.mgrname, a.type, a.rdate, a.date,
    sum(a.shares_adj*b.price)/1000000 as Assets
       label="Total Portfolio Assets, x$1m",
from Holding as a, crsp_m as b
where a.permno=b.permno and a.rdate=b.qdate
group by a.mgrno, a.rdate
order by mgrno, rdate;
quit;
```

## 4. Quarterly mutual funds holdings data (mf_holding.csv & mf_asset.csv)
```
%let mfldate = 01JAN2012;

/* Get report and vintage dates from Thomson-Reuters Mutual Fund Holdings */
/* Exclude Non-Equity Funds from Holdings data that is reported as of Fiscal Quarter End */
/* First, Keep First Vintage with Holdings Data for Each RDATE-FUNDNO */
proc sql;
  create table First_Vint
  as select distinct intnx("month",rdate,0,"E") as rdate format YYMMDD10., fdate, fundno
  from tfn2.s12type1
  where ("&begdate"d <= rdate <="&enddate"d and ioc not in (1,5,6,7))
  group by fundno, intnx("month",rdate,0,"E")
  having fdate=min(fdate) and max(rdate)=rdate
  order by fundno, rdate desc; 
quit; 

data First_Vint; set First_Vint;  
  by fundno descending rdate;	format nrdate YYMMDD10.;
  nrdate = lag(rdate); if first.fundno then nrdate = intnx("month",rdate,6,"E");
run; 

proc sort data=First_Vint nodupkey; by fundno fdate; run; 

/* Add WFICN portfolio identifiers from MFLINKS */
proc sql;
  create table First_Vint
  as select b.wficn as wficn_old, a.*
  from First_Vint as a left join  mfl.mflink2 as b
  on  a.fundno=b.fundno and a.fdate = b.fdate;
quit; 

data First_Vint; 
  length WFICN 6.;
  set First_Vint;
  by fundno fdate; retain wficn;
  if first.fundno or not missing(wficn_old) then wficn=wficn_old;
  if missing(wficn_old) and fdate<"&mfldate"d then wficn=.;
  drop wficn_old;
run;

proc sort data=First_Vint nodupkey; by wficn rdate; run; 

/* Extract Holdings by Merging TR-MF S12type1 and S12type3 Sets */
proc sql;
  create view MF_Holdings  /* Add Holdings Data */
    as select a.rdate,a.nrdate,a.fdate,a.wficn,a.fundno,b.cusip,b.shares
    from First_Vint as a, tfn2.s12type3 as b
    where a.fdate=b.fdate and a.fundno=b.fundno; 
  create view MF_Holdings1 /* Map Holdings CUSIP to CRSP Unique Identifier PERMNO */
    as select a.rdate,a.nrdate,a.fdate,a.wficn,a.fundno,b.permno,a.shares
    from MF_Holdings as a, (select distinct ncusip, permno from crsp.msenames
                         where not missing(ncusip)) as b
    where a.cusip=b.ncusip; 
  create table MF_Holdings2 /* Adjust Shares on Vintage Dates */
    as select a.rdate,a.nrdate, a.wficn, a.fundno, a.permno,
        a.shares*b.cfacshr as shares_adj label = "Adjusted Shares Held"
    from MF_Holdings1 as a, crsp_m as b
    where a.permno=b.permno and a.fdate=b.qdate;  
quit;

proc sort data=MF_Holdings2 nodupkey thread; by wficn permno rdate; where shares_adj>0; run;

proc sql;
  create table MF_Assets
  as select distinct a.wficn, a.rdate, 
      sum(a.shares_adj*b.price)/1000000 as Assets
        label="Total Portfolio Assets, x$1m"
  from MF_Holdings2 as a, crsp_m as b
  where a.permno=b.permno and a.rdate=b.qdate
  group by a.wficn, a.rdate
  order by wficn, rdate;
quit;
```

## 5. Domestic equity mutual funds return & tna data (mf_returns.csv)
```
%let md = (year(date)*12 + month(date) - 1925*12 - 11);

/* CRSP policies to exclude: note that it is available through 1990 only */
%let crsp_policy_to_exclude='C & I','Bal','Bonds','Pfd','B & P','GS','MM','TFM';

/* lipper classes to include */
%let lipper_class=('EIEI','G','LCCE','LCGE','LCVE','MCCE','MCGE','MCVE',
                           'MLCE','MLGE','MLVE','SCCE','SCGE','SCVE');

%let lipper_obj_cd=('CA','EI','G','GI','MC','MR','SG');

/* strategic insight classes to include */
%let si_obj_cd=('AGG','GMC','GRI','GRO','ING','SCG');

/* weisenberger classes to include */
%let wbrger_obj_cd=('G','GCI','IEQ','LTG','MCG','SCG');

/* retain latest style for every crsp_fundno */
data fund_style; set crsp.fund_style;
      format lipper_classX $4.lipper_obj_cdX $3. si_obj_cdX $3. wbrger_obj_cdX $5.;
      by crsp_fundno;
      %macro temp (var=);
            retain &var.X;
            if first.crsp_fundno then &var.X = &var;
            else if &var ne '' then &var.X = &var;
      %mend temp;
      %temp(var=lipper_class);
      %temp(var=lipper_obj_cd);
      %temp(var=si_obj_cd);
      %temp(var=wbrger_obj_cd);
run;

 
data latest_style (rename=(lipper_classX=lipper_class lipper_obj_cdX=lipper_obj_cd
      si_obj_cdX=si_obj_cd wbrger_obj_cdX=wbrger_obj_cd)); set fund_style;
      drop lipper_class lipper_obj_cd si_obj_cd wbrger_obj_cd;
      by crsp_fundno;
      if last.crsp_fundno;
      keep crsp_fundno lipper_classX lipper_obj_cdX si_obj_cdX wbrger_obj_cdX crsp_obj_cd;
run;

/* keep most recent name of every share class */
data names; set crsp.fund_names (where=(fund_name ne ''));
      by crsp_fundno;
      if last.crsp_fundno;
      keep ncusip ticker crsp_fundno fund_name index_fund_flag mgmt_cd mgmt_name first_offer_dt;
run;

 
/* domestic equity mutual funds: this borrows from return_gap.sas sample program on
   WRDS, http://goo.gl/DB1iz (you have to be a WRDS member to access that link) */

data funds; set fund_style (keep=crsp_fundno si_obj_cd wbrger_obj_cd policy lipper_class
      lipper_obj_cd where=(policy not in (&crsp_policy_to_exclude)));
      long_way = 1;
      dummy = 0;
      if lipper_class in &lipper_class or lipper_obj_cd in (&lipper_obj_cd) then output;
      else if missing(lipper_class)=1 and missing(lipper_obj_cd)=1 and
            si_obj_cd in &si_obj_cd then output;
      else if missing(lipper_class)=1 and missing(lipper_obj_cd)=1 and missing(si_obj_cd)=1
            and wbrger_obj_cd in &wbrger_obj_cd then output;

      /* note that these are unnecessary: for all but seven funds, crsp_fundnos are
         identified as equity funds without these statements, and none of the seven
         funds that have all dummy = 1 (and no dummy = 0) cases is a domestic
         diversified equity funds, so actually that merging to get avrcs is not needed */

      /* else if missing(lipper_class)=1 and missing(lipper_obj_cd)=1 and missing(si_obj_cd)=1
            and missing(wbrger_obj_cd)=1 and policy='CS' then do;
            dummy = 1; output;end;

      else if missing(lipper_class)=1 and missing(lipper_obj_cd)=1 and missing(si_obj_cd)=1
            and missing(wbrger_obj_cd)=1 and missing(policy)=1 and 80<=avrcs<=105 then do;
            dummy = 1; output; end; */

      keep crsp_fundno long_way;

run;

proc sort data = funds nodupkey; by crsp_fundno; run;

/* shorter way of identifying objectives using crsp_obj_cd variable */
data funds2; set crsp.fund_style;
      if substr(crsp_obj_cd,1,1) = 'E'; /* equity */
      if substr(crsp_obj_cd,2,1) = 'D'; /* domestic */
      if substr(crsp_obj_cd,3,1) in ('C','Y'); /* cap-based or style */
      if substr(crsp_obj_cd,3,2) not in ('YH','YS'); /* exclude hedged, short */
      if si_obj_cd ne 'OPI'; /* exclude option income */
      short_way = 1;
      keep crsp_fundno short_way;
run;

proc sort data = funds2 nodupkey; by crsp_fundno; run;

/* there are only 85 funds that don't overlap between long way and short way: all
   of them are in the CRSP sample (the one obtained using crsp_obj_cd) but not in mine:
   all of them appear to have changed style to international or are otherwise international */

data funds3; merge funds funds2 latest_style names;
      by crsp_fundno;
      if long_way = 1 or short_way = 1;
      if long_way = short_way then delete;
run;

/* look for funds that have flip-flopped their style */
proc sql;
      create table funds4 as select a.*, b.crsp_obj_cd
      from funds as a left join fund_style (keep=crsp_fundno crsp_obj_cd begdt) as b
      on a.crsp_fundno = b.crsp_fundno
      order by a.crsp_fundno, begdt;
quit;

data funds4; set funds4;
      if substr(crsp_obj_cd,1,3) not in ('EDC','EDY') then flipper = 1;
      else flipper = 0;
run;

proc sql;
      create table funds5 as select crsp_fundno, max(flipper) as flipper from funds4
      group by crsp_fundno;
quit;

/* identify index and target date funds and drop them from the sample */
data funds6; merge funds5 (in=in1) names;
      by crsp_fundno;
      if in1;
run;

proc sql;
      create table funds7 (drop=index_fund_flag) as select * from funds6 group by crsp_fundno
      having count(index_fund_flag) = 0;
quit;

data funds8; set funds7;
      format namex $140.;
      namex = lowcase(fund_name);
      if max(index(namex,'index'), index(namex,'s&p')) = 0;
      if max(index(namex,'idx'), index(namex,'dfa'), index(namex,'program')) = 0;
      if max(indexw(namex,'etf'),index(namex,'exchange traded'),
            index(namex,'exchange-traded')) = 0;
      if max(index(namex,'target'),index(namex,'2005'),index(namex,'2005'),
            index(namex,'2010'),index(namex,'2015'),index(namex,'2020'),index(namex,'2025'),
            index(namex,'2030'),index(namex,'2035'),index(namex,'2040'),
            index(namex,'2045'),index(namex,'2050'),index(namex,'2055')) = 0;
      drop namex ;
run;

/* group into funds using mflinks and exclude flippers: drops sample slightly. there are
   some unusual observations where some share classes change style while other share
   classes in the same fund do not (e.g., wficn=100001) */
proc sql;
      create table equity_funds (where=(wficn ne .)) as
      select b.wficn, a.fund_name, a.ncusip, a.ticker, a.crsp_fundno, a.mgmt_cd, a.mgmt_name, a.first_offer_dt from funds8 as a left join mfl.mflink1 as b
      on a.crsp_fundno=b.crsp_fundno
      group by b.wficn having max(flipper) = 0 order by crsp_fundno;
quit;

/* delete unneeded datasets */
proc datasets nowarn nolist nodetails;
      delete funds funds2-funds8 fund_style latest_style;
run; quit;

 
/********************************************************************************************

      PART 3: create a dataset that for every equity fund WFICN and month dummy MD contains size

      TNA, net return RET, and gross return RRET. Some code below borrows from return_gap.sas sample

      program on WRDS, http://goo.gl/DB1iz (you have to be a WRDS member to access that link)

********************************************************************************************/

data returns1; merge equity_funds (in=in1) crsp.monthly_tna_ret_nav (drop=mnav);
      by crsp_fundno;
      if in1 = 1;
      if wficn ne .;
      if mtna < 0 then mtna = .;
      mtna = mtna + 0;
      mret = mret + 0;
      /* retain last available tna */
      retain tna;
      if first.crsp_fundno or mtna ne . then tna = mtna;
      rename caldt=date;
      drop mtna;
run;

/* get expense ratio data required to get gross returns */
proc sql;
      create table returns2 as select a.*, b.exp_ratio, b.turn_ratio
      from returns1 (rename=(tna=mtna)) as a left join crsp.fund_fees as b
      on a.crsp_fundno=b.crsp_fundno and date between b.begdt and b.enddt
      order by crsp_fundno, date;
quit;

/* compute gross returns */
data returns2; set returns2;
      by crsp_fundno date;
      if exp_ratio = -99 then exp_ratio=.;
      weight = lag(mtna);
      if first.crsp_fundno then weight = mtna;
      rret=sum(mret,exp_ratio/12);
run;

/* aggregate multiple share classes */
proc sort data = returns2; by wficn date; run;

data multiclass1 oneclass; set returns2;
      by wficn date;
      if first.date=0 or last.date=0 then output multiclass1;
      else output oneclass;
run;

proc sql;
      create table multiclass2
      as select wficn, date, mgmt_name, sum(mret*weight)/sum(weight) as mret,
      sum(mtna) as mtna, sum(rret*weight)/sum(weight) as rret,
      sum(exp_ratio*weight)/sum(weight) as exp_ratio, sum(turn_ratio*weight)/sum(weight) as turn_ratio,
      Min(first_offer_dt) as first_offer_dt, FIRST(mgmt_cd) as mgmt_cd
      from multiclass1 group by wficn, date;
quit;

data returns; set oneclass (drop=weight) multiclass2;
      md = &md;
      rename mret=ret mtna=tna;
run;

proc sort data = returns nodupkey; by wficn md; run;

proc datasets nowarn nolist nodetails;
      delete returns1 ;
quit;

```

## 6. Fama and French's Book to Market ratio
[Open Source Asset Pricing](https://www.openassetpricing.com/data/)

## 7. NYSE ME breakpoints data (ME_breakpoints.csv)
[Kenneth R. French's Data Library](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html)