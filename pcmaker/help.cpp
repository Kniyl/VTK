// help.cpp : implementation file
//

#include "stdafx.h"
#include "pcmaker.h"
#include "help.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// help dialog


help::help(CWnd* pParent /*=NULL*/)
	: CDialog(help::IDD, pParent)
{
	//{{AFX_DATA_INIT(help)
	//}}AFX_DATA_INIT
	m_helptext = 
    "Help for pcmaker.exe\r\n"
    "====================\r\n\r\n"
    "pcmaker.exe is a program that generates Makefiles "
    "for vtk.  This is the first step in the process of "
    "compiling vtk on a PC. There are a few pieces of "
    "information that are required in order to build vtk. "
    "The first is the directory location of your vtk "
    "source code. This should not have any spaces in it. "
    "C:\\vtk2.0 is an example of what someone might type. "
    "Next you need to tell pcmaker where to place the "
    "resulting object files and libraries. C:\\vtkbin is "
    "one example. The directory you enter does not have "
    "to exist, pcmaker will ask you if you want it created.\r\n\r\n"
    "Next you'll need to provide the directory where your "
    "C++ compiler is installed. Some examples are:\r\n"
    "  C:\\msdev\r\n"
    "  C:\\program Files\\DevStudio\\vc\r\n"
    "  D:\\msvc40\r\n"
    "This path can have a space in it.\r\n\r\n"
    "If you plan to build Java libraries then you need to "
    "specify the path to your copy of the Java JDK. If "
    "you do not plan on building java libraries then you "
    "should leave this blank. Otherwise pcmaker will "
    "complain that it cannot find the required header files."
    ;
}


void help::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(help)
	DDX_Text(pDX, IDC_EDIT1, m_helptext);
	DDV_MaxChars(pDX, m_helptext, 10000);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(help, CDialog)
	//{{AFX_MSG_MAP(help)
		// NOTE: the ClassWizard will add message map macros here
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// help message handlers
