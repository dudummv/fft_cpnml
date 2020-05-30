
(*Funções para dividir os vetores*)

fun split2 (x1::x2::xs) =
    (case split2 xs of
          (ys, zs) => (x1::ys, x2::zs))
  | split2 rest = (rest, [])

fun split3 (x1::x2::x3::xs) =
    (case split3 xs of
          (ys, zs, ts) => (x1::ys, x2::zs, x3::ts))
       
  | split3 (rest::rest2) = ([rest],rest2 , [])
  | split3 rest = (rest,[],[])

fun split5 (x1::x2::x3::x4::x5::xs) =
    (case split5 xs of
          (ys, zs, ts, ws, rs) => (x1::ys, x2::zs, x3::ts, x4::ws, x5::rs))
       
  | split5 (rest::rest2::rest3::rest4) = ([rest], [rest2] , [rest3], rest4, [])
  | split5 (rest::rest2::rest3) = ([rest], [rest2] , rest3, [], [])
  | split5 (rest::rest2) = ([rest],rest2 ,[], [], [])
  | split5 rest = (rest,[] ,[], [], [])

(*Converte de (r) => (r,0) *)

fun realToImag x =
  List.map(fn k => (k,0.0)) x

(*Produto entre complexos*)
fun complexProd ((r:real,i:real),(rr:real,ii:real)) =
  (r*rr-i*ii,r*ii+i*rr)
(*Soma entre complexos*)
fun complexSum ((r:real,i:real),(rr:real,ii:real)) = 
  (r+rr,i+ii)

(*Conjugado complexo*)
fun conjugate(r,i) =
  (r,~1.0*i)

(*Função auxiliar dos calculos*)
fun join(x,y,w) =
  if List.null(x) then
    []
  else 
    complexSum(List.hd(x),complexProd(List.hd(y),List.hd(w)))::join(List.tl(x),List.tl(y),List.tl(w))

(*Produto entre listas de complexos*)
fun complexProdList(x,y) =
  if List.null(x) then
    []
  else
    complexProd(List.hd(x),List.hd(y))::complexProdList(List.tl(x),List.tl(y))

fun calcWphase(k,n,const,mlt,phase) =
  (Math.cos(const*Math.pi*real(k)/(mlt*n) + phase), ~1.0*Math.sin(const*Math.pi*real(k)/(mlt*n) + phase))


(*fft radix 2, utilizada no algoritimo de bluestein*)

fun fft2 c = 
  if List.length(c) <> 1 then
    let
      val (l1,l2) = split2(c)
      val x1 = fft2(l1)
      val x2 = fft2(l2)
      
      val n = List.length(x1)
      
      val w1 = List.tabulate(n,(fn k => calcWphase(k,real(n),2.0,2.0,0.0)))
     
      val r1 = join(x1,x2,w1)
      val r2 = join(x1,x2,(List.map(fn x => complexProd(x,(~1.0,0.0))) w1))
    in
      r1@r2
    end  
  else
    c
(*ifft radix 2, utilizada no algoritimo de bluestein*)
fun ifft2 c =
  let
    val x = List.map conjugate c
    val n = real(List.length(x))
    val y = fft2(x)
  in
    List.map(fn (r,i) => (r/n,~1.0*i/n)) y
  end

(*MRFFT*)
fun fft c =
  (*Radix 5*)
  if List.length(c) mod 5 = 0 then
    let
      val (l1,l2,l3,l4,l5) = split5(c)
      val x1 = fft(l1)
      val x2 = fft(l2)
      val x3 = fft(l3)
      val x4 = fft(l4)
      val x5 = fft(l5) 

      val n = List.length(x1)
      (*Constantes*)
      val t1 = 2.0*Math.pi/5.0
      val t2 = 4.0*Math.pi/5.0
      val t3 = 6.0*Math.pi/5.0
      val t4 = 8.0*Math.pi/5.0

      val w1 = List.tabulate(n,(fn k => calcWphase(k,real(n),2.0,5.0,0.0)))
      val w2 = List.tabulate(n,(fn k => calcWphase(k,real(n),4.0,5.0,0.0)))
      val w3 = List.tabulate(n,(fn k => calcWphase(k,real(n),6.0,5.0,0.0)))
      val w4 = List.tabulate(n,(fn k => calcWphase(k,real(n),8.0,5.0,0.0)))


      val w5 = List.tabulate(n,(fn k => calcWphase(k,real(n),2.0,5.0,t1)))
      val w6 = List.tabulate(n,(fn k => calcWphase(k,real(n),4.0,5.0,t2)))
      val w7 = List.tabulate(n,(fn k => calcWphase(k,real(n),6.0,5.0,t3)))
      val w8 = List.tabulate(n,(fn k => calcWphase(k,real(n),8.0,5.0,t4)))

      val w9 = List.tabulate(n,(fn k => calcWphase(k,real(n),2.0,5.0,t2)))
      val w10 = List.tabulate(n,(fn k => calcWphase(k,real(n),4.0,5.0,t4)))
      val w11 = List.tabulate(n,(fn k => calcWphase(k,real(n),6.0,5.0,t1)))
      val w12 = List.tabulate(n,(fn k => calcWphase(k,real(n),8.0,5.0,t3)))

      val w13 = List.tabulate(n,(fn k => calcWphase(k,real(n),2.0,5.0,t3)))
      val w14 = List.tabulate(n,(fn k => calcWphase(k,real(n),4.0,5.0,t1)))
      val w15 = List.tabulate(n,(fn k => calcWphase(k,real(n),6.0,5.0,t4)))
      val w16 = List.tabulate(n,(fn k => calcWphase(k,real(n),8.0,5.0,t2)))

      val w17 = List.tabulate(n,(fn k => calcWphase(k,real(n),2.0,5.0,t4)))
      val w18 = List.tabulate(n,(fn k => calcWphase(k,real(n),4.0,5.0,t3)))
      val w19 = List.tabulate(n,(fn k => calcWphase(k,real(n),6.0,5.0,t2)))
      val w20 = List.tabulate(n,(fn k => calcWphase(k,real(n),8.0,5.0,t1)))

      
      val r1 = join(join(join(join(x1,x2,w1),x3,w2),x4,w3),x5,w4)
      val r2 = join(join(join(join(x1,x2,w5),x3,w6),x4,w7),x5,w8)
      val r3 = join(join(join(join(x1,x2,w9),x3,w10),x4,w11),x5,w12)
      val r4 = join(join(join(join(x1,x2,w13),x3,w14),x4,w15),x5,w16)
      val r5 = join(join(join(join(x1,x2,w17),x3,w18),x4,w19),x5,w20)
    in
      r1@r2@r3@r4@r5
    end 
  (*Radix 3*)
  else if List.length(c) mod 3 = 0 then
    let
      val (l1,l2,l3) = split3(c)
      val x1 = fft(l1)
      val x2 = fft(l2)
      val x3 = fft(l3)
      
      val n = List.length(x1)
      (*Constantes*)
      val t1 = 2.0*Math.pi/3.0
      val t2 = 4.0*Math.pi/3.0

      val w1 = List.tabulate(n,(fn k => calcWphase(k,real(n),2.0,3.0,0.0)))
      val w2 = List.tabulate(n,(fn k => calcWphase(k,real(n),4.0,3.0,0.0)))
      val w3 = List.tabulate(n,(fn k => calcWphase(k,real(n),2.0,3.0,t1)))
      val w4 = List.tabulate(n,(fn k => calcWphase(k,real(n),4.0,3.0,t2)))
      val w5 = List.tabulate(n,(fn k => calcWphase(k,real(n),2.0,3.0,t2)))
      val w6 = List.tabulate(n,(fn k => calcWphase(k,real(n),4.0,3.0,t1)))
      val r1 = join(join(x1,x2,w1),x3,w2)
      val r2 = join(join(x1,x2,w3),x3,w4)
      val r3 = join(join(x1,x2,w5),x3,w6)   
    in
      r1@r2@r3
    end  
  (*Radix 2*)
  else if List.length(c) mod 2 = 0 then
    let
      val (l1,l2) = split2(c)
      val x1 = fft(l1)
      val x2 = fft(l2)
      
      val n = List.length(x1)
      
      val w1 = List.tabulate(n,(fn k => calcWphase(k,real(n),2.0,2.0,0.0)))
     
      val r1 = join(x1,x2,w1)
      val r2 = join(x1,x2,(List.map(fn x => complexProd(x,(~1.0,0.0))) w1))
    in
      r1@r2
    end
  (*Caso Base*)   
  else if List.length(c) = 1 then
    c
  (*Bluestein*)  
  else
    let
      (*Proxima potencia para radix*)
      val nextPow = Real.ceil(Math.log10(2.0*real(List.length(c))-1.0)/Math.log10(2.0))
      val newN = Math.pow(2.0,real(nextPow))
      val n = List.length(c)
     (*Sequencias usadas na convolução*)
      val bN1 = List.tabulate(n,(fn k => Math.pow(real(k),2.0)/2.0))
      val bN2 = List.tabulate(n-1,(fn k => Math.pow(real(k-n+1),2.0)/2.0))@bN1
         
      val aa1 = List.map(fn k => calcWphase(1,real(n),2.0*k,1.0,0.0)) bN1
      val aa2 = List.map(fn k => calcWphase(1,real(n),2.0*k,1.0,0.0)) bN2
      

      val naa = List.length(aa2)
      (*Sinal com zero padding*)
      val newX = complexProdList(c,aa1)@List.tabulate (Real.round(newN)-n, (fn _ => (0.0,0.0)))
      (*Sequencia com zero padding*)
      val y = (List.map(fn k => conjugate(k)) aa2)@List.tabulate (Real.round(newN)-naa, (fn _ => (0.0,0.0)))
      (*FFT de radix2*)
      val fy = fft2(y)
      val fx = fft2(newX)
      (*Convolução no tempo*)
      val res = complexProdList(fx,fy)
      val res2 = ifft2(res)
      (*Retirando apenas a parte significativa do sinal*)
      val g1 = List.take(List.drop(res2,n-1),n)
    in
      (*Multiplicação final*)
     complexProdList(g1,aa1)
    end

(*IFFT usando a tecnica da FFT direta no conjugado complexo*)
fun ifft c =
  let
    val x = List.map conjugate c
    val n = real(List.length(x))
    val y = fft(x)
  in
    List.map(fn (r,i) => (r/n,~1.0*i/n)) y
  end

  

