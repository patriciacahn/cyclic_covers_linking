class Link_Component:
    """
    Attributes
    ----------
    num_components: number of components
    """
    def __init__(self, num_components: int, signs=[], overstrands=[]):
        self.num_components=num_components
        self.signs=signs
        self.overstrands=overstrands
        
        

class Link_Diagram:
    """
    Attributes
    ----------
    num_components: number of components
    branch_component: which component plays the role of the branch curve
    """
    def __init__(self, num_components: int, branch_component=0):
        self.num_components=num_components
        self.branch_component=0
        self.component_list=[]
        
    def add_component(self, component : Link_Component):
        self.component_list.append(component)
        
    
class Cyclic_Cover:
    """
    Attributes
    ----------
    num_components: number of components
    branch_component: which component plays the role of the branch curve
    """
    def __init__(self, degree : int, diagram : Link_Diagram):
        self.branch_component=diagram.branch_component
        self.degree=degree
        self.diagram=diagram
        self.C=CyclicPermutationGroup(degree)
        self.cycle=self.C(tuple([i+1 for i in range(degree)]))
    
    def three_cell_list(self,component_number):
        three_cell_list=[]
        three_cell_list.append(self.C.identity())
        for i in range(len(self.diagram.component_list[component_number].signs)):
            if self.diagram.component_list[component_number].signs[i]==1 and self.diagram.component_list[component_number].overstrands[i][0]==self.branch_component:
                three_cell_list.append(self.cycle*three_cell_list[i])
            elif self.diagram.component_list[component_number].signs[i]==-1 and self.diagram.component_list[component_number].overstrands[i][0]==self.branch_component:
                three_cell_list.append(self.cycle^-1*three_cell_list[i])
            else:
                three_cell_list.append(three_cell_list[i])
        return three_cell_list
    
    def two_chain(self,component_number,lift_number):
        #component_number is which pseudo-branch curve P we are using
        #lift_number is the lift P^q of P that we are bounding
        
        #2-chain Sigma^q= \sum B_i^q + \sum x_i^j(A_i^j), superscripts j taken mod p
        #i is the crossing number
        #A_i^1... A_i^p are lifts of vertical 2-cell A_i
        
        #n = number of arcs of K. i \in {0,...,n-1}
        #j \in {1,...,p}
        #n*p variables x_i^j
        #p equations per crossing, plus another n equations given by \sum_j x_i^j=0 for fixed i \in {0,...,n-1}
        #n*p+n equations
        #Store system of equations in (pn+n) x (pn) matrix (not augmented)
        
        n=len(self.diagram.component_list[self.branch_component].signs)
        p=self.degree
        
       
        
        M=Matrix(n*p+n,n*p)
        #columns x_0^1, x_0^2,\dots, x_0^p, x_1^1,x_1^2,\dots,x_1^p,\dots
        #p equations for each crossing
        #rows 0 to p-1 are equations for crossing 0, etc
        #final n rows are the equations \sum_j x_i^j=0 for fixed i \in {0,...,n-1}
       
        #vector of constant terms
        v=vector([0 for i in range(n*p+n)])
        
        #M|v is the augmented matrix
        
        
        def column_index(i,j):
            return (p*i+j-1)%(n*p)
        def row_index(i,j):
            return p*i+j-1
        
        knot_overstrands=self.diagram.component_list[self.branch_component].overstrands
        knot_signs=self.diagram.component_list[self.branch_component].signs
        
        #loop through strands of the knot 
        for i in range(n):
            #identify the crossing type
            
            #overstrand=knot
            if knot_overstrands[i][0]==self.branch_component:
                
                for j in range(1,p+1):
                    #x_i^j-x_{i+1}^j-e(i)(x_{o(i)}^s-x_{o(i)}^{t})=0, t=s+1 %p  between {1,...,p}
                    #w(o(i))
                    w_over=self.three_cell_list(self.branch_component)[knot_overstrands[i][1]]                        
                    #w(i)
                    w_under=self.three_cell_list(self.branch_component)[i]
                    
                    #superscript of cell above; depends on crossing sign
                    if knot_signs[i]==1:
                        s=w_over.inverse()(w_under(j))
                    elif knot_signs[i]==-1:
                        if w_under(j)!=1:
                            s=w_over.inverse()(w_under(j)-1)
                        else:
                            s=w_over.inverse()(p)
                    
                    t=s+1
                    if t==p+1:
                        t=1
                    
                    #coefficient of x_i^j is 1
                    M[row_index(i,j),column_index(i,j)]+=1
                    #coefficient of x_{i+1}^j is -1
                    M[row_index(i,j),column_index(i+1,j)]+=-1
                    #coefficient of x_{o(i)}^s is -e(i)
                    M[row_index(i,j),column_index(knot_overstrands[i][1],s)]+=-knot_signs[i]
                    #coefficient of x_{o(i)}^t is e(i)
                    M[row_index(i,j),column_index(knot_overstrands[i][1],t)]+=knot_signs[i]
                    
                
            #overstrand=pb curve whose qth lift we are bounding  
            elif knot_overstrands[i][0]==component_number:
                
                for j in range(1,p+1):
                    #x_i^j-x_{i+1}^j=e(i) x_i^(j-1)-x_{i+1}^(j-1)=-e(i) if s_j=q
                    #x_i^j-x_{i+1}^j=0 else
                    #w(o(i)) 
                    w_over=self.three_cell_list(component_number)[knot_overstrands[i][1]]
                    #w(i)
                    w_under=self.three_cell_list(self.branch_component)[i]
                    s_j=w_over.inverse()(w_under(j))#which lift of the pb curve is at page (j,j+1) of the roladex
                    
                    
                    #coefficient of x_i^j is 1
                    M[row_index(i,j),column_index(i,j)]+=1
                    #coefficient of x_{i+1}^j is -1
                    M[row_index(i,j),column_index(i+1,j)]+=-1
                    #constant term
                    if s_j==lift_number:
                        v[row_index(i,j)]+=knot_signs[i]
                        if j!=1:
                            v[row_index(i,j-1)]+=-knot_signs[i]
                        else:
                            v[row_index(i,p)]+=-knot_signs[i]
                        
                        
                    
             
            #overstrand=another pb curve 
            else:
                for j in range(1,p+1):
                    #x_i^j-x_{i+1}^j=0 
                    #coefficient of x_i^j is 1
                    M[row_index(i,j),column_index(i,j)]+=1
                    #coefficient of x_{i+1}^j is -1
                    M[row_index(i,j),column_index(i+1,j)]+=-1
                    
            
            #for each knot arc i, have the equation x_i^1+...+x_i^p=0
            #these go in the last n rows of the np+n x np matrix
            #note last n rows of v (constant terms) are already 0
            
            for j in range(1,p+1):
                M[n*p+i,column_index(i,j)]+=1
                
            
            
        
        #Uncomment to display equations
        #for i in range(n*p+n):
        #    print(M.row(i))
        #    for j in range(n*p):
        #        if M[i,j]!=0:
        #            print("+",M[i,j],"*x_",j // p,"^",(j% p)+1)
        #    print("=",v.column()[i])
        
        
        return M.solve_right(v)
    
    def linking_number(self,component_number_1,lift_number_1,component_number_2,lift_number_2):
        #Compute linking number of component_number_1, lift_number_1 with component_number_2,lift_number_2 by intersecting 
        #component_number_2,lift_number_2 with chain bounding component_number_1, lift_number_1
        
        pb_1_chain=self.two_chain(component_number_1,lift_number_1)
        
        pb_2_overstrands=self.diagram.component_list[component_number_2].overstrands
        pb_2_signs=self.diagram.component_list[component_number_2].signs
        
        linking_number=0
        
        for i in range(len(pb_2_overstrands)):
            arc=pb_2_overstrands[i]
            crossing_type=arc[0]
            over_arc_number=arc[1]
            
            
            if crossing_type==component_number_1:
                
                    #w(o(i))
                    w_over=self.three_cell_list(component_number_1)[over_arc_number]
                    #w(i)
                    w_under=self.three_cell_list(component_number_2)[i]
                    s=w_over.inverse()(w_under(lift_number_2))
                    
                    if s==lift_number_1:
                        linking_number+=pb_2_signs[i]
                        #print(pb_2_signs[i])
                
            elif crossing_type==self.branch_component:
                    #w(o(i))
                    w_over=self.three_cell_list(self.branch_component)[over_arc_number]
                    #print("w_over")
                    #print(w_over)
                    #w(i)
                    w_under=self.three_cell_list(component_number_2)[i]
                    #print("w_under")
                    #print(w_under)
                    #superscript of cell above; depends on crossing sign
                    if pb_2_signs[i]==1:
                        s=w_over.inverse()(w_under(lift_number_2))
                    elif pb_2_signs[i]==-1:
                        if w_under(lift_number_2)!=1:
                            s=w_over.inverse()(w_under(lift_number_2)-1)
                        else:
                            s=w_over.inverse()(p)
                    
                    #value of x_o(i)^s * sign
                    linking_number+=pb_2_signs[i]*pb_1_chain[self.degree*over_arc_number+s-1]
                    #print("s",s)
                    #print("over arc",over_arc_number)
                    #print("index",self.degree*over_arc_number+s-1)
                    #print(pb_2_signs[i]*pb_1_chain[self.degree*over_arc_number+s-1])
                    
                
            #elif crossing_type==component_number_2:
        return linking_number
            


        
        
        
    
    
    