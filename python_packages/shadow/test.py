#===================== test ======================
from shadow.plot import *


def fit_example():
    r'''

    def my_func(x, a, b, c):
        return (a*x**2 + b*x + c)*np.sin(x)**2

    label = '${}x^2+{}x+{}$'

    xs = np.arange(-4,4,0.1)

    ys = my_func(xs,1,1,1) + np.random.random(len(xs))

    line_plot(xs,ys,'s',label='Scatter')

    #=========fit==========
    _, text1 = fit(xs,ys,xs,':',func=my_func,precision=1,label=label)

    #=========labels=======
    title('Title')
    xlabel('x-axis')
    ylabel('y-axis')
    legend_on()

    '''
    set_markersize(8)

    def my_func(x, a, b, c):
        return (a*x**2 + b*x + c)*np.sin(x)**2

    label = '${}x^2+{}x+{}$'

    xs = np.arange(-4,4,0.1)

    ys = my_func(xs,1,1,1) + np.random.random(len(xs))

    line_plot(xs,ys,'s',label='Scatter')

    #=========fit==========
    _, text1 = fit(xs,ys,xs,':',func=my_func,precision=1,label=label)

    #=========labels=======
    title('Title')
    xlabel('x-axis')
    ylabel('y-axis')
    legend_on()


def energy_landscape_example():
    fig,[ax] = panel(1,1,l_p=0,r_p=0,b_p=0,t_p=0,dpi=100)
    rectangle(-1000,1000,-1000,1000,fc='g',ec='m',lw=4)
    plt.axis('off')

    xlim([-1000,1000])
    ylim([-1000,1000])

    from scipy.stats import norm

    np.random.seed(584111)

    x = np.arange(-1000,1000,1)
    y = 50*np.sin(0.02*x)+50*np.sin(0.05*x)-300*np.sin(0.007*x)+500
    for i in range(1,11):
        y -= np.random.random()*2000*norm.pdf((-1000+300*i+x)/20)
    y[0],y[-1] = -1000,-1000
    ax.fill(x,y,c='m',alpha=0.2)
    ax.plot(x,y,c='m',alpha=0.5)

    pos = [-670]
    for i in pos:
        xn = np.argmax(i<x)
        ax.arrow(x[xn]+30+100,y[xn]+30,-50,0,width=15,fc='grey',ec='grey')
        ax.arrow(x[xn]-30-100,y[xn]+30,50,0,width=15,fc='grey',ec='grey')
        ellipse([x[xn],y[xn]+30],20,alpha=1,fc='b')

    pos = [320]
    for i in pos:
        xn = np.argmax(i<x)
        ax.arrow(x[xn]+30,y[xn]+30,50,0,width=15,fc='grey',ec='grey')
        ax.arrow(x[xn]-30,y[xn]+30,-50,0,width=15,fc='grey',ec='grey')
        ellipse([x[xn],y[xn]+30],20,alpha=1,fc='b')

    text(400,y[1400]-100,r'\textbf{Global minima}',ha='center')
    text(250+100,y[1250]-300,r'\textbf{Local minima}',ha='right')
    ax.arrow(225,y[1225]-200,0,130,width=15,fc='k',)
    arrow(pos=[[-300,y[1250]-300],[x[330],y[330]]],curve=-0.5,)



def pdf_example():
    r'''
        fig,[ax] = panel(1,1,l_p=0,r_p=0,b_p=0,t_p=0,dpi=100)
        plt.axis('off')


        # add circles
        ellipse([0,0],200,fc='none',linewidth=8,ec='r',)
        ellipse([0,0],350,fc='none',linewidth=10,ec='g')
        ellipse([0,0],400,fc='none',linewidth=10,ec='b')


        #make a hexagon (dashed line)
        obj1 = polygon([200,0],200,sides=6,fc='none',ls=':',ec='k',lw='2')

        data = obj1.get_path().vertices
        obj = [obj1]

        # make 2 more hexagon (dashed line)
        obj2 = polygon([-100,data[1][1]],200,sides=6,fc='none',ls=':',ec='k',lw='2')
        obj3 = polygon([-100,-data[1][1]],200,sides=6,fc='none',ls=':',ec='k',lw='2')

        # make 3 hexagon (hidden) to use their vertices as centers for atoms
        obj.append(polygon([-400,0],200,sides=6,fc='none'))
        obj.append(polygon([-data[0][0]+data[1][0],3*data[1][1]],200,sides=6,fc='none'))
        obj.append(polygon([-data[0][0]+data[1][0],-3*data[1][1]],200,sides=6,fc='none'))

        # add atoms
        for data in obj:
            for i in data.get_path().vertices[:-1]:
                if sum(i**2)<450**2:
                    polygon(i+0*(-0.5+np.random.random()),20,sides=20,fc='k',alpha=1,ec='k',lw=2)

        ax.arrow(0,550,0,550,width=1,head_width=20,fc='k')
        ax.arrow(-50,600,850,0,width=1,head_width=20,fc='k')

        ax.arrow(200,-400,0,350,width=1,head_width=20,fc='k')
        ax.arrow(200,-400,400,0,width=1,head_width=1,fc='k')

        ax.arrow(350,-250,0,200,width=1,head_width=20,fc='k')
        ax.arrow(350,-250,250,0,width=1,head_width=1,fc='k')

        ax.arrow(400,-100,0,50,width=1,head_width=20,fc='k')
        ax.arrow(400,-100,200,0,width=1,head_width=1,fc='k')

        set_font({'size':15})

        text(650,-400,r'\textbf{1st Coordination Shell}')
        text(650,-250,r'\textbf{2nd Coordination Shell}')
        text(650,-100,r'\textbf{3rd Coordination Shell}')

        set_font({'size':20})

        line_plot([180,180],[0,600],'r:',lw=1,alpha=0.4)
        line_plot([220,220],[0,600],'r:',lw=1,alpha=0.4)


        line_plot([330,330],[0,600],'g:',lw=1,alpha=0.4)
        line_plot([370,370],[0,600],'g:',lw=1,alpha=0.4)

        line_plot([380,380],[0,600],'b:',lw=1,alpha=0.4)
        line_plot([420,420],[0,600],'b:',lw=1,alpha=0.4)

        line_plot([180,190,200,210,220],[600,950,1000,950,600],'r',lw=1)
        line_plot([330,340,350,360,370],[600,1000,1200,1000,600],'g',lw=1)
        line_plot([380,390,400,410,420],[600,750,800,750,600],'b',lw=1)

        text(-100,900,'$g(r)$',rotation=90)
        text(300,540,'$r$')

        xlim([-600,1200])
        ylim([-600,1200])


    '''
    fig,[ax] = panel(1,1,l_p=0,r_p=0,b_p=0,t_p=0,dpi=100)
    plt.axis('off')


    # add circles
    ellipse([0,0],200,fc='none',linewidth=8,ec='r',)
    ellipse([0,0],350,fc='none',linewidth=10,ec='g')
    ellipse([0,0],400,fc='none',linewidth=10,ec='b')


    #make a hexagon (dashed line)
    obj1 = polygon([200,0],200,sides=6,fc='none',ls=':',ec='k',lw='2')

    data = obj1.get_path().vertices
    obj = [obj1]

    # make 2 more hexagon (dashed line)
    obj2 = polygon([-100,data[1][1]],200,sides=6,fc='none',ls=':',ec='k',lw='2')
    obj3 = polygon([-100,-data[1][1]],200,sides=6,fc='none',ls=':',ec='k',lw='2')

    # make 3 hexagon (hidden) to use their vertices as centers for atoms
    obj.append(polygon([-400,0],200,sides=6,fc='none'))
    obj.append(polygon([-data[0][0]+data[1][0],3*data[1][1]],200,sides=6,fc='none'))
    obj.append(polygon([-data[0][0]+data[1][0],-3*data[1][1]],200,sides=6,fc='none'))

    # add atoms
    for data in obj:
        for i in data.get_path().vertices[:-1]:
            if sum(i**2)<450**2:
                polygon(i+0*(-0.5+np.random.random()),20,sides=20,fc='k',alpha=1,ec='k',lw=2)

    ax.arrow(0,550,0,550,width=1,head_width=20,fc='k')
    ax.arrow(-50,600,850,0,width=1,head_width=20,fc='k')

    ax.arrow(200,-400,0,350,width=1,head_width=20,fc='k')
    ax.arrow(200,-400,400,0,width=1,head_width=1,fc='k')

    ax.arrow(350,-250,0,200,width=1,head_width=20,fc='k')
    ax.arrow(350,-250,250,0,width=1,head_width=1,fc='k')

    ax.arrow(400,-100,0,50,width=1,head_width=20,fc='k')
    ax.arrow(400,-100,200,0,width=1,head_width=1,fc='k')

    set_font({'size':15})

    text(650,-400,r'\textbf{1st Coordination Shell}')
    text(650,-250,r'\textbf{2nd Coordination Shell}')
    text(650,-100,r'\textbf{3rd Coordination Shell}')

    set_font({'size':20})

    line_plot([180,180],[0,600],'r:',lw=1,alpha=0.4)
    line_plot([220,220],[0,600],'r:',lw=1,alpha=0.4)


    line_plot([330,330],[0,600],'g:',lw=1,alpha=0.4)
    line_plot([370,370],[0,600],'g:',lw=1,alpha=0.4)

    line_plot([380,380],[0,600],'b:',lw=1,alpha=0.4)
    line_plot([420,420],[0,600],'b:',lw=1,alpha=0.4)

    line_plot([180,190,200,210,220],[600,950,1000,950,600],'r',lw=1)
    line_plot([330,340,350,360,370],[600,1000,1200,1000,600],'g',lw=1)
    line_plot([380,390,400,410,420],[600,750,800,750,600],'b',lw=1)

    text(-100,900,'$g(r)$',rotation=90)
    text(300,540,'$r$')

    xlim([-600,1200])
    ylim([-600,1200])
