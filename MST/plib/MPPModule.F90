!*************************************************************************!
! NAME: MPPModule.F90                                                     !
!                                                                         !
! HISTORY ---------                                                       !
!                                                                         !
!     VERSION: 2.5b                                                       !
!     AUTHOR: Yang Wang (07/09/2018)                                      !
!     ADDED FEATURE: Added one function for allowing the MPI processes    !
!                    to perform non-blocking send package                 !
!                    nbsendPackage: non-block sending a package           !
!            Note: nbrecvPackage is not implemented, because of a major   !
!                  techinical difficulty due to that the package size is  !
!                  not known beforehand and needs to call MPI_iprobe to   !
!                  find out. However, even though MPI_iprobe is non-      !
!                  blocking, it doesn't return a MPI_request handler for  !
!                  MPI_wait to use.                                       !
!                                                                         !
!     VERSION: 2.5a                                                       !
!     AUTHOR: Yang Wang (01/31/2017)                                      !
!     ADDED FEATURE: Added two functions for allowing the MPI processes   !
!                    to take action in turns:                             !
!                    startRoundTurn:  start a round robin tuin count      !
!                    finishMyTurn:    finish my turn in the round robin   !
!                                                                         !
!     VERSION: 2.1a                                                       !
!     AUTHOR: Yang Wang (04/05/2013)                                      !
!     ADDED FEATURE: Define MPP_BUFFER_MEM, a guess of MPI buffer size.   !
!                    In allreduce calls, MPI_IN_PLACE is not used if the  !
!                    array size if greater than the buffer size.          !
!                                                                         !
!     VERSION: 2.0k                                                       !
!     AUTHOR: Yang Wang (09/14/2011)                                      !
!     ADDED FEATURE: add the following function:                          !
!                    getMaxWaits: get the maximum number of nonblocking   !
!                                 messages before waitMessage is called   !
!                    setMaxWaits: set the maximum number of nonblocking   !
!                                 messages before waitMessage is called   !
!                                                                         !
!     VERSION: 2.0j                                                       !
!     AUTHOR: Yang Wang (08/31/2011)                                      !
!     ADDED FEATURE: add the following two functions:                     !
!                    isThereMessage: returns True/False if there is a     !
!                                    message arrived.                     !
!                                                                         !
!     VERSION: 2.0i                                                       !
!     AUTHOR: Yang Wang (06/09/2011)                                      !
!     ADDED FEATURE: Use MPI_IN_PLACE for the send buffer in MPI_Allgather!
!                    routines to comply with the latest MPI implementation!
!                    requirement.                                         !
!                                                                         !
!     VERSION: 2.0h                                                       !
!     AUTHOR: Aurelian Rusanu & Yang Wang (09/13/2010)                    !
!     ADDED FEATURE: Got rid off PVM code, which is considered obsolete,  !
!                    added NumPEs_comm and MyPE_comm, and modified        !
!                    setCommunicator, resetCommunicator, and changed      !
!                    NumPEs and MyPE in communication routines to         !
!                    NumPEs_comm and MyPE_comm, respectively.             !
!                                                                         !
!     VERSION: 2.0g                                                       !
!     AUTHOR: Markus Eisenbach (12/11/2006)                               !
!     ADDED FEATURE: provide a MPI communicator at initialization.        !
!                 * original call to MPPinit() will set the communicator  !
!                   to MPI_COMM_WORLD.                                    !
!                 * MPPinit(comm) will cause the MPP module to use the    !
!                   specified communicator.                               !
!                                                                         !
!     VERSION: 2.0f                                                       !
!     ADDED FEATURE: chang LongIntKind for message types to IntKind       !
!                                                                         !
!     VERSION: 2.0e                                                       !
!     AUTHIOR: Yang Wang (04/28/2005)                                     !
!     ADDED FEATURE: overload openLocalMemory to allow multi-dimension    !
!                    local array to be read by remote processes.          !
!                                                                         !
!     VERSION: 2.0d                                                       !
!     AUTHIOR: Yang Wang (04/20/2005)                                     !
!     ADDED FEATURE: add a function returning MPI_COMM_WORLD:             !
!                    getCommunicator                                      !
!                                                                         !
!     VERSION: 2.0c                                                       !
!     AUTHIOR: Yang Wang (05/04/2004)                                     !
!     ADDED FEATURE: LAM MPI library complains MPI_REAL8 not defined.     !
!                    I changed MPI_REAL8 to MPI_DOUBLE_PRECISION.         !
!                    The current LAM MPI does not support MPI_win_lock and!
!                    MPI_win_unlock. I add MPI_win_fence subroutine call  !
!                    so that if No_MPI_LOCK is true, MPI_win_fence will be!
!                    called instead of MPI_win_lock/MPI_win_unlock.       !
!                    Therefore, option -DNo_MPI_LOCK needs to be added to !
!                    the preprocessor if the MPI library does not support !
!                    MPI_win_lock/MPI_win_unlock.                         !
!                  * I have to add syncLocalMemory as a functionality     !
!                    that allows MPI_win_fence to be called. This routine !
!                    will do nothing if the MPI library is not LAM.       !
!                    supported.                                           !
!                  * I also noticed that MPI_get does not work for        !
!                    complex data type in LAM, so I use real data type    !
!                    instead with twice as much the array size. The change!
!                    is applied to readRemoteMemory_c. I suspect that     ! 
!                    MPI_put will not work properly for complex data      !
!                    type in LAM neither, so I also use real data type    !
!                    instead with twice as much the array size. The change!
!                    is applied to writeRemoteMemory_c.                   ! 
!                                                                         !
!     VERSION: 2.0b                                                       !
!     AUTHIOR: Yang Wang (04/29/2004)                                     !
!     ADDED FEATURE: PGI compiler complains argument mismatch between the !
!                    strcopy and its calling routine: pack_str1 and       !
!                    pack_str2. The following routines are modified to    !
!                    fix this problem:                                    !
!                    pack_str1, pack_str2, strcopy                        !
!                                                                         !
!                                                                         !
!     VERSION: 2.0a                                                       !
!     AUTHIOR: Yang Wang (07/13/2003)                                     !
!     ADDED FEATURE: Fixed bugs in the following routines for one-sided   !
!                    communications:                                      !
!                    readRemoteMemory  : in MPI_get, CHARACTER_BYTES,     !
!                                        INTEGER_BYTES, REAL_BYTES, and   !
!                                        COMPLEX_BYTES are changed to     !
!                                        MPI_CHARACTER, MPI_INTEGER,      !
!                                        MPI_DOUBLE_PRECISION, and        !
!                                        MPI_DOUBLE_COMPLEX respectively. !
!                    writeRemoteMemory : in MPI_put, CHARACTER_BYTES,     !
!                                        INTEGER_BYTES, REAL_BYTES, and   !
!                                        COMPLEX_BYTES are changed to     !
!                                        MPI_CHARACTER, MPI_INTEGER,      !
!                                        MPI_DOUBLE_PRECISION, and        !
!                                        MPI_DOUBLE_COMPLEX respectively. !
!                                                                         !
!     VERSION: 2.0                                                        !
!     AUTHIOR: Yang Wang (10/14/2002)                                     !
!     ADDED FEATURE: The following routines for one-sided communications  !
!                    openLocalMemory   : make a local memory space to     !
!                                        allow remote processes to access !
!                    readRemoteMemory  : copy over the contents of the    !
!                                        memory space from a remote       !
!                                        processor                        !
!                    writeRemoteMemory : write to the memory space of a   !
!                                        remote processor                 !
!                    closeLocalMemory  : make the local memory space,     !
!                                        which is previously opened for   !
!                                        remote access, unaccessible      !
!                                                                         !
!     VERSION: 1.3c                                                       !
!     AUTHIOR: Yang Wang (04/10/2003)                                     !
!     MODIFIED: MaxWaits                                                  !
!               MaxWaits is changed from:                                 !
!                  integer (kind=IntKind), parameter :: MaxWaits=5        !
!               to:                                                       !
!                  integer (kind=IntKind), parameter :: MaxWaits=200      !
!                                                                         !
!     VERSION: 1.3b                                                       !
!     AUTHIOR: Yang Wang (04/10/2003)                                     !
!     MODIFIED: sendPackage0, sendPackage1, recvPackage                   !
!          Add the following lines to the routines to avoid possible      !
!          breaking down in the case of NumPEs = 1:                       !
!                                                                         !
!           else if (.not.isParentAlive() .and. NumPEs.eq.1) then         !
!              return                                                     !
!                                                                         !
!                                                                         !
!     VERSION: 1.3a                                                       !
!     AUTHIOR: Yang Wang (05/18/2002)                                     !
!           In pack_str2, the following bug is fixed:                     !
!              size=n*m*l                                                 !
!           is changed to                                                 !
!              size=n*m*la                                                !
!           Eliminated flush(6) from the code to avoid error message      !
!           on CRAY-T3E                                                   !
!           Fixed bug in calling MPI_allgather where the last 2 arguments:!
!           status,info are replace with 1 argument: info.                !
!                                                                         !
!     VERSION: 1.3                                                        !
!     AUTHIOR: Yang Wang (04/24/2002)                                     !
!     ADDED FEATURE: initParent                                           !
!                    startChildTasks : must be called by parent           !
!                    isParent  : true if calling process is parent        !
!                    isChild   : true if calling process is child         !
!                    getParent : return Parent PE                         !
!                    getNumPEs : return NumPEs excluding parent           !
!                    getMyPE   : return MyPE (=0, 1, ...)                 !
!                                if I am parent, MyPE = NumPEs            !
!                    getHostPEs: return PE-id array on a host machine     !
!                                This can only be called by parent        !
!                    killParent: kill parent process                      !
!                    endParent : need to call killParent first            !
!                    isParentAlive: true, if Parent is not killed         !
!     MODIFIED: initMPP                                                   !
!               changed original private function getNumPEs to            !
!               subroutine calNumPEs                                      !
!                                                                         !
!     VERSION: 1.2c                                                       !
!     AUTHIOR: Yang Wang (04/09/2002)                                     !
!     MODIFIED: Changed MaxWaits from 100 to 5                            !
!     ADDED FEATURE: printWaitInfo                                        !
!                    getNumWaits                                          !
!                                                                         !
!     VERSION: 1.2b                                                       !
!     AUTHIOR: Yang Wang (10/17/2001)                                     !
!     MODIFIED: GlobalCollect routines for MPI using MPI_allgather        !
!                                                                         !
!     VERSION: 1.2a                                                       !
!     AUTHIOR: Yang Wang (05/04/2000)                                     !
!     BUGS FIXED:  1. in unpackMessage, added bpos=0 to initialize bpos.  !
!                     this is necessary for CRAY MPP systems              !
!                  2. in glb_collect_*, use MPI_isend instead of MPI_send !
!                     to avoid possible job hung due to multiple message  !
!                     send                                                !
!                  3. implement PVM version of nbsendMessage, waitMessage,!
!                     and nbrecvMessage (05/10/2000)                      !
!                                                                         !
!     VERSION: 1.2                                                        !
!     AUTHIOR: Yang Wang (02/24/2000)                                     !
!     ADDED FEATURE: using LongIntKind for message types                  !
!                    nbsendMessage                                        !
!                    nbrecvMessage                                        !
!                    waitMessage                                          !
!                                                                         !
!     VERSION: 1.1                                                        !
!     AUTHIOR: Yang Wang (11/04/1999)                                     !
!     ADDED FEATURE: packMessage                                          !
!                    unpackMessage                                        !
!                    sendPackage                                          !
!                    recvPackage                                          !
!                                                                         !
!     VERSION: 1.0 - beta                                                 !
!     AUTHIOR: Yang Wang (10/16/1999)                                     !
!     FEATURE: The original design                                        !
!                                                                         !
!                                                                         !
! DISCRIPTION -----                                                       !
!     A Fortran 90 module.                                                !
!     It contains message passing interfaces and supports the following   !
!     subroutine calls. The actual message passing libaray specific to    !
!     the system is hidden from the user.                                 !
!     * MyPE:                                                             !
!       an integer variable that equals to the current node number        !
!     * NumPEs:                                                           !
!       an integer variable that equals to the number of nodes being used !
!       Note: parent is not counted                                       !
!     * AllPEs:                                                           !
!       an integer parameter used for sending a message to all PEs        !
!     * AnyPE:                                                            !
!       an integer parameter used for receiving a message from any PE     !
!     * initMPP():                                                        !
!       initialize message passing module and must be called first        !
!     * endMPP():                                                         !
!       deallocate and terminate the message passing module, and is       !
!       recommended to be called before the end                           !
!     * getSource():                                                      !
!       returns the source node that sended the message, which was just   !
!       received by the current node.                                     !
!     * syncAllPEs():                                                     !
!       Synchronize all the nodes                                         !
!                                                                         !
!     * sendMessage(message,msgtype,target_node):                         !
!       send message to the target node                                   !
!       message: a value of integer, real, complex, or character type     !
!       msgtyp: message id number                                         !
!       target_node: the node that the message is sent to.                !
!                    If = -1, the message is sent to all the nodes.       !
!     * sendMessage(message,size1,msgtype,target_node):                   !
!       send 1-d array message to the target node                         !
!       message: 1-d array of integer, real, complex, or character type   !
!       size1: the size of the 1-d array                                  !
!       msgtyp: message id number                                         !
!       target_node: the node that the message is sent to                 !
!                    If = -1, the message is sent to all the nodes.       !
!     * sendMessage(message,size1,size2,msgtype,target_node):             !
!       send 2-d array message to the target node                         !
!       message: 2-d array of integer, real, complex, or character type   !
!       size1: the size of the first dimension of the array               !
!       size2: the size of the second dimension of the array              !
!       msgtyp: message id number                                         !
!       target_node: the node that the message is sent to                 !
!                    If = -1, the message is sent to all the nodes.       !
!                                                                         !
!     * msendMessage(message,msgtype,target_nodes,node_counts)            !
!       send message to multiple target nodes                             !
!       message: a value of integer, real, complex, or character type     !
!       msgtyp: message id number                                         !
!       target_nodes: the array of nodes that the message is sent to.     !
!       node_count: the number of target nodes.                           !
!     * msendMessage(message,size1,msgtype,target_nodes,node_counts)      !
!       send 1-d array message to multiple target nodes                   !
!       message: 1-d array of integer, real, complex, or character type   !
!       size1: the size of the 1-d array                                  !
!       msgtyp: message id number                                         !
!       target_nodes: the array of nodes that the message is sent to.     !
!       node_count: the number of target nodes.                           !
!     * msendMessage(message,size1,size2,msgtype,target_nodes,node_counts)!  
!       send 2-d array message to multiple target nodes                   !
!       message: 2-d array of integer, real, complex, or character type   !
!       size1: the size of the first dimension of the array               !
!       size2: the size of the second dimension of the array              !
!       msgtyp: message id number                                         !
!       target_nodes: the array of nodes that the message is sent to.     !
!       node_count: the number of target nodes.                           !
!                                                                         !
!     * recvMessage(message,msgtype,source_node):                         !
!       receive message from the source node                              !
!       message: a value of integer, real, complex, or character type     !
!       msgtyp: message id number                                         !
!       source_node: the node that receives message                       !
!                    If = -1, the message is received from any node       !
!     * recvMessage(message,size1,msgtype,source_node):                   !
!       receive 1-d array message from the source node                    !
!       message: 1-d array of integer, real, complex, or character type   !
!       size1: the size of the 1-d array                                  !
!       msgtyp: message id number                                         !
!       source_node: the node that receives the message                   !
!                    If = -1, the message is received from any node       !
!     * recvMessage(message,size1,size2,msgtype,source_node):             !
!       receive 2-d array message from the source node                    !
!       message: 2-d array of integer, real, complex, or character type   !
!       size1: the size of the first dimension of the array               !
!       size2: the size of the second dimension of the array              !
!       msgtyp: message id number                                         !
!       source_node: the node that receives the message                   !
!                    If = -1, the message is received from any node       !
!                                                                         !
!     * bcastMessage(message,root_node):                                  !
!       broadcast the message to all the nodes by the root node.          !
!       must be called by all the nodes in the same time.                 !
!       if MyPE is not root_node, message is a returned value.            !
!       message: a value of integer, real, complex, or character type     !
!       root_node: the node that sends the message                        !
!     * bcastMessage(message,size1,root_node):                            !
!       broadcast the 1-d array message to all the nodes by the root node.!
!       must be called by all the nodes in the same time.                 !
!       if MyPE is not root_node, message is a returned value.            !
!       message: 1-d array of integer, real, complex, or character type   !
!       size1: the size of the 1-d array                                  !
!       root_node: the node that sends the message                        !
!     * bcastMessage(message,size1,size2,root_node):                      !
!       broadcast the 2-d array message to all the nodes by the root node.!
!       must be called by all the nodes in the same time.                 !
!       if MyPE is not root_node, message is a returned value.            !
!       message: 2-d array of integer, real, complex, or character type   !
!       size1: the size of the first dimension of the array               !
!       size2: the size of the second dimension of the array              !
!       root_node: the node that sends the message                        !
!                                                                         !
!     * GlobalMin(message):                                               !
!     * GlobalMax(message):                                               !
!     * GlobalSum(message):                                               !
!       Calculate the minimum, maximum, or summation of the message across!
!       all the nodes. The final result is stored in message. It should be!
!       called by all the nodes.                                          !
!       message: a value of integer, real, complex, or character type     !
!     * GlobalMin(message,size1):                                         !
!     * GlobalMax(message,size1):                                         !
!     * GlobalSum(message,size1):                                         !
!       Calculate the minimum, maximum, or summation of the 1-d message   !
!       array across all the nodes. The final result is stored in the     !
!       message array. It should be called by all the nodes.              !
!       size1: the size of the 1-d array                                  !
!     * GlobalMin(message,size1,size2):                                   !
!     * GlobalMax(message,size1,size2):                                   !
!     * GlobalSum(message,size1,size2):                                   !
!       Calculate the minimum, maximum, or summation of the 2-d message   !
!       array across all the nodes. The final result is stored in the     !
!       message array. It should be called by all the nodes.              !
!       size1: the size of the first dimension of the array               !
!       size2: the size of the second dimension of the array              !
!     * GlobalSum(message, size1, size2, size3)                           !
!       Calculate the summation of the 3-d message array across all nodes.!
!       The final result is stored in the message array. It should be     !
!       called by all nodes.                                              !
!       size1: the size of the first dimension of the array               !
!       size2: the size of the second dimension of the array              !
!       size3: the size of the third dimension of the array               !
!                                                                         !
!     * GlobalCollect(message):                                           !
!       Collect the value message(MyPE) from all the node and store in    !
!       the array message(1:NumPEs). It should be called by all the nodes.!
!     * GlobalCollect(message, size1):                                    !
!       Collect the 1-d array of message(1:size1,MyPE) from all the node  !
!       and store in the array message(1:size1,1:NumPEs). It should be    !
!       called by all the nodes.                                          !
!       size1: the size of the 1-d array                                  !
!     * GlobalCollect(message, size1, size2):                             !
!       Collect the 2-d array of message(1:size1,1:size2,MyPE) from all   !
!       the node and store in the array message(1:size1,1:size2,1:NumPEs).!
!       It should be called by all the nodes.                             !
!       size1: the size of the first dimension of the array               !
!       size2: the size of the second dimension of the array              !
!                                                                         !
!     * packMessage(message):                                             !
!       Pack the (integer, string, real, or complex) message to be sent   !
!     * packMessage(message, size1)                                       !
!       Pack the (integer, string, real, or complex) 1-d array message to !
!       be send                                                           !
!       size1: the size of the 1-d array                                  !
!     * packMessage(message, size1, size2):                               !
!       Pack the (integer, string, real, or complex) 2-d array message to !
!       be send                                                           !
!       size1: the size of the first dimension of the array               !
!       size2: the size of the second dimension of the array              !
!                                                                         !
!     * unpackMessage(message):                                           !
!       Un-pack the (integer, string, real, or complex) message           !
!     * unpackMessage(message, size1)                                     !
!       Un-pack the (integer, string, real, or complex) 1-d array message !
!       size1: the size of the 1-d array                                  !
!     * unpackMessage(message, size1, size2):                             !
!       Un-pack the (integer, string, real, or complex) 2-d array message !
!       size1: the size of the first dimension of the array               !
!       size2: the size of the second dimension of the array              !
!                                                                         !
!     * sendPackage(msgtype,target_node):                                 !
!       send the packed messages to the target node                       !
!       if target_node = AllPEs, the messages are sent to all PEs         !
!     * sendPackage(msgtype,target_nodes,num_of_nodes):                   !
!       send the packed messages to the target nodes specified by the     !
!       array target_nodes of size num_of_nodes                           !
!     * recvPackage(msgtype,source_node):                                 !
!       receive the packed messages from the source node                  !
!       if source_node = AnyPE, the node can be any PE                    !
!                                                                         !
!     * nbsendMessage(message,msgtype,target_node):                       !
!     * nbsendMessage(message,size,msgtype,target_node):                  !
!     * nbsendMessage(message,size1,size2,msgtype,target_node):           !
!                                                                         !
!     * nbrecvMessage(message,msgtype,source_node):                       !
!     * nbrecvMessage(message,size,msgtype,source_node):                  !
!     * nbrecvMessage(message,size1,size2,msgtype,source_node):           !
!                                                                         !
!     * waitMessage(MsgID):                                               !
!       wait for the message whose id is MsgID                            !
!     * printWaitInfo():                                                  !
!       print out information about the messages that are waiting to be   !
!       sent or received                                                  !
!     * getNumWaits():                                                    !
!       return an integer number that is the number of messages on the    !
!       waiting list                                                      !
!                                                                         !
!     * openLocalMemory(space,size,accessID)                              !
!       make a local memory space to allow remote process to access       !
!       space(:) :: input, 1-d array for open to remote access            !
!       size     :: input, integer, size of 1-d array space               !
!       accessID :: output, integer, access ID                            !
!       Note: This is a synchronized call.                                !
!     * readRemoteMemory(localArray,n,remoteProc,accessID,startIndex):    !
!       copy over the contents of the memory space of size n starting     !
!       from index startIndex from remote processor remoteProc to local   !
!       array localArray.                                                 !
!       localArray(:) :: output, 1-d array containing contents copied     !
!                        from the memory space on a remote processor      !
!       n             :: input, integer, size of the contents to be       !
!                        copied over from remote processor                !
!       remoteProc    :: input, integer, remote processor ID              !
!       accessID      :: input, integer, access ID of the remote memory   !
!       startIndex    :: input, integer, the starting index of the        !
!                        remote memory for reading, assuming that the     !
!                        first index of the remote array is 1.            !
!     * writeRemoteMemory(localArray,n,remoteProc,accessID,startIndex):   !
!       write the contents of the local array of size n to remote memory  !
!       space starting from index starting from index startIndex on       !
!       processor remoteProc.                                             !
!       localArray(:) :: input, 1-d array containing contents to be       !
!                        written to the memory space on a remote processor!
!       n             :: input, integer, size of the contents to be       !
!                        written to the memory space on a remote processor!
!       remoteProc    :: input, integer, remote processor ID              !
!       accessID      :: input, integer, access ID of the remote memory   !
!       startIndex    :: input, integer, the starting index of the        !
!                        remote memory for writing, assuming that the     !
!                        first index of the remote array is 1.            !
!     * closeLocalMemory(AccessID):                                       !
!       make the local memory space, which is previously opened to allow  !
!       remote access, unaccessible                                       !
!       accessID :: input, integer, access ID                             !
!       Note: This is a synchronized call.                                !
!     * syncLocalMemory(AccessID):                                        !
!       synchronize the local memory space, which has been opened for     !
!       remote access. This should be called after calls of read/write    !
!       remote meomory. This helps to update the local memory.            !
!       accessID :: input, integer, access ID                             !
!       Note: This is a synchronized call.                                !
!                                                                         !
!     * setComunicator(comm,p,n):                                         !
!       set the communicator to comm, MyPE to p, and NumPEs to n, so that !
!       all the communication functions will work for the new universe    !
!       defined by comm.                                                  !
!     * resetCommunicator():                                              !
!       set the communicator back to the default one (MPI_COMM_WORLD).    !
!*************************************************************************!
      module MPPModule
         use KindParamModule, only : IntKind, RealKind, CmplxKind, LongIntKind
         use ErrorHandlerModule, only : StopHandler
         use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
         implicit none
#ifdef MPI
         include 'mpif.h'
         integer status(MPI_STATUS_SIZE)
#else
#define MPI_ADDRESS_KIND IntKind
#define MPI_COMM_WORLD 0
#endif
!
         integer (kind=IntKind), parameter :: INTEGER8 = 8
         integer (kind=IntKind), parameter :: INTEGER4 = 4
!
         interface initMPP
           module procedure initMPP_default_comm, initMPP_comm
         end interface
!
         interface sendMessage
            module procedure send_int0, send_real0, send_cmplx0
            module procedure send_int1, send_real1, send_cmplx1
            module procedure send_int2, send_real2, send_cmplx2
            module procedure send_str0, send_str1,  send_str2
         end interface
!
         interface msendMessage
            module procedure msend_int0, msend_real0, msend_cmplx0
            module procedure msend_int1, msend_real1, msend_cmplx1
            module procedure msend_int2, msend_real2, msend_cmplx2
            module procedure msend_str0, msend_str1,  msend_str2
         end interface
!
         interface bcastMessage
            module procedure bcast_int0, bcast_real0, bcast_cmplx0
            module procedure bcast_int1, bcast_real1, bcast_cmplx1
            module procedure bcast_int2, bcast_real2, bcast_cmplx2
            module procedure bcast_str0, bcast_str1,  bcast_str2
         end interface
!
         interface recvMessage
            module procedure recv_int0, recv_real0, recv_cmplx0
            module procedure recv_int1, recv_real1, recv_cmplx1
            module procedure recv_int2, recv_real2, recv_cmplx2
            module procedure recv_str0, recv_str1,  recv_str2
         end interface
!
         interface GlobalMax
            module procedure glb_max_int0, glb_max_real0, glb_max_cmplx0
            module procedure glb_max_int1, glb_max_real1, glb_max_cmplx1
            module procedure glb_max_int2, glb_max_real2, glb_max_cmplx2
         end interface
!
         interface GlobalMin
            module procedure glb_min_int0, glb_min_real0, glb_min_cmplx0
            module procedure glb_min_int1, glb_min_real1, glb_min_cmplx1
            module procedure glb_min_int2, glb_min_real2, glb_min_cmplx2
         end interface
!
         interface GlobalSum
            module procedure glb_sum_int0, glb_sum_real0, glb_sum_cmplx0
            module procedure glb_sum_int1, glb_sum_real1, glb_sum_cmplx1
            module procedure glb_sum_int2, glb_sum_real2, glb_sum_cmplx2
            module procedure glb_sum_int3, glb_sum_real3, glb_sum_cmplx3
         end interface
!
         interface nbsendMessage
            module procedure nbsend_int0, nbsend_real0, nbsend_cmplx0
            module procedure nbsend_int1, nbsend_real1, nbsend_cmplx1
            module procedure nbsend_int2, nbsend_real2, nbsend_cmplx2
            module procedure nbsend_str0, nbsend_str1,  nbsend_str2
         end interface
!
         interface nbrecvMessage
            module procedure nbrecv_int0, nbrecv_real0, nbrecv_cmplx0
            module procedure nbrecv_int1, nbrecv_real1, nbrecv_cmplx1
            module procedure nbrecv_int2, nbrecv_real2, nbrecv_cmplx2
            module procedure nbrecv_str0, nbrecv_str1,  nbrecv_str2
         end interface
!
         interface GlobalCollect
            module procedure glb_collect_int0, glb_collect_real0,   &
     &                       glb_collect_cmplx0, glb_collect_str0
            module procedure glb_collect_int1, glb_collect_real1,   &
     &                       glb_collect_cmplx1, glb_collect_str1
            module procedure glb_collect_int2, glb_collect_real2,   &
     &                       glb_collect_cmplx2, glb_collect_str2
         end interface
!
         interface packMessage
            module procedure pack_int0, pack_real0, pack_cmplx0,   &
     &                       pack_str0, pack_longint0
            module procedure pack_int1, pack_real1, pack_cmplx1,   &
     &                       pack_str1, pack_longint1
            module procedure pack_int2, pack_real2, pack_cmplx2,   &
     &                       pack_str2, pack_longint2
         end interface
!
         interface unpackMessage
            module procedure unpack_int0, unpack_real0, unpack_cmplx0,   &
     &                       unpack_str0, unpack_longint0
            module procedure unpack_int1, unpack_real1, unpack_cmplx1,   &
     &                       unpack_str1, unpack_longint1
            module procedure unpack_int2, unpack_real2, unpack_cmplx2,   &
     &                       unpack_str2, unpack_longint2
         end interface
!
         interface sendPackage
            module procedure sendPackage0, sendPackage1
         end interface
!
         interface openLocalMemory
            module procedure openLocalMemory_s, openLocalMemory_i,   &
     &                       openLocalMemory_r, openLocalMemory_c
            module procedure openLocalMemory_s1, openLocalMemory_i1,   &
     &                       openLocalMemory_r1, openLocalMemory_c1
         end interface
!
         interface readRemoteMemory
            module procedure readRemoteMemory_s, readRemoteMemory_i, &
     &                       readRemoteMemory_r, readRemoteMemory_c
         end interface
!
         interface writeRemoteMemory
            module procedure writeRemoteMemory_s, writeRemoteMemory_i, &
     &                       writeRemoteMemory_r, writeRemoteMemory_c
         end interface
!
         public :: initMPP
         public :: endMPP
         public :: getSource
         public :: getMyPE
         public :: getNumPEs
         public :: getHostPEs
!
         public :: bcastMessage
         public :: recvMessage
         public :: sendMessage
         public :: msendMessage
         public :: packMessage
         public :: unpackMessage
         public :: recvPackage
         public :: sendPackage
         public :: nbsendPackage
         public :: syncAllPEs
         public :: GlobalCollect
         public :: GlobalSum
         public :: GlobalMax
         public :: GlobalMin
         public :: nbsendMessage
         public :: nbrecvMessage
!
         public :: initParent
         public :: isParent
         public :: isChild
         public :: isParentAlive
         public :: startChildTasks
         public :: getParent
         public :: endParent
!
         public :: waitMessage
         public :: printWaitInfo
         public :: getNumWaits
!
         public :: openLocalMemory
         public :: readRemoteMemory
         public :: writeRemoteMemory
         public :: closeLocalMemory
         public :: syncLocalMemory
!
         public :: getCommunicator
         public :: setCommunicator
         public :: resetCommunicator
!
         public :: isThereMessage
         public :: getMaxWaits
         public :: setMaxWaits
!
         public :: startRoundTurn
         public :: finishMyTurn

         integer (kind=IntKind), public :: MPP_BUFFER_MEM = 1024000 ! A conservative guess
!
         integer (kind=IntKind), public :: NumPEs = 1
         integer (kind=IntKind), public :: MyPE  = 0
         integer (kind=IntKind), parameter, public :: AllPEs = -1
#ifdef MPI
         integer (kind=IntKind), parameter, public :: AnyPE = MPI_ANY_SOURCE
#else
         integer (kind=IntKind), parameter, public :: AnyPE = -1
#endif
!
         private
         integer (kind=IntKind), parameter :: MaxIOPEs=32
!
         integer (kind=IntKind) :: MyInstance
         integer :: NumChildTasks
!
         integer, allocatable :: TaskID(:)
         integer :: ParentPE = -2
         integer :: ParentTaskID = -1
         integer :: NumHosts = 0
!
         integer :: info
         integer (kind=IntKind) :: All = AllPEs
!
!        the communicator to use for MPI default will be MPI_COMM_WORLD
         integer (kind=IntKind) :: communicator
!
         integer (kind=IntKind) :: i, j, k, slen
         integer (kind=IntKind) :: sourceNode
         integer (kind=IntKind) :: bufferID
!
         integer (kind=IntKind), parameter :: MaxWaits_param = 50000
         integer (kind=IntKind) :: MaxWaits = MaxWaits_param
         integer (kind=IntKind) :: NumWaits
         integer (kind=IntKind) :: NumSends, NumRecvs
         integer (kind=IntKind), allocatable :: MsgID(:)
         integer (kind=IntKind), allocatable :: MsgBuffer(:)
         integer (kind=IntKind), allocatable :: MsgSRC(:)
         character (len=4), allocatable :: WaitType(:)
!
!        integer (kind=IntKind) :: ai
!        integer (kind=IntKind), allocatable :: ai1(:)
!        integer (kind=IntKind), allocatable :: ai2(:,:)
!
!        real (kind=RealKind) :: ar
!        real (kind=RealKind), allocatable :: ar1(:)
!        real (kind=RealKind), allocatable :: ar2(:,:)
!
!        complex (kind=CmplxKind) :: ac
!        complex (kind=CmplxKind), allocatable :: ac1(:)
!        complex (kind=CmplxKind), allocatable :: ac2(:,:)
!
!        character, allocatable :: as(:)
!        character, allocatable :: as1(:,:)
!        character, allocatable :: as2(:,:,:)
!
!=====================================================================
!        Pack/unpack related parameters, types, and variables
!=====================================================================
!
         logical :: received
         logical :: packing
!
         integer, parameter :: MaxNumPacks = 5000
!
         integer :: NumPacks
!
#ifdef MPI
         type BufferPointers
            integer size
            integer DataType
            character (len=1), pointer :: sBuf(:)
            integer (kind=IntKind), pointer :: iBuf(:)
            real (kind=RealKind), pointer :: rBuf(:)
            complex (kind=CmplxKind), pointer :: cBuf(:)
         end type BufferPointers
!
         character (len=1), allocatable :: bspace(:)
!
         type (BufferPointers), save :: PackedBuffer(MaxNumPacks)
#endif
!
         integer :: bsize, bpos
         integer, public :: INTEGER_BYTES
         integer, public :: REAL_BYTES
         integer, public :: COMPLEX_BYTES
         integer, public :: CHARACTER_BYTES
         integer, public :: ADDRESS_BYTES
!
         logical :: Parent = .false.
         logical :: Child = .false.
         logical :: ChildTasksStarted = .false.
         logical :: ParentKilled = .false.
!
         type HostType
             integer :: NumPEs
             character (len=100) :: address
             character (len=100) :: executable
             integer, pointer :: PEs(:)
         end type HostType
         type (HostType), allocatable :: Hosts(:)
         integer :: MyHostIndex = 0
!
         integer (kind=IntKind) :: MyPE_comm = 0
         integer (kind=IntKind) :: NumPEs_comm = 1
!
         integer (kind=IntKind) :: robin_message_type = 20170131
         integer (kind=IntKind) :: first_robin, robin_wait
!
      contains
!=====================================================================
!=====================================================================
!=====================================================================
      subroutine initParent()
         implicit none
         Parent = .true.
         NumPEs = 0
         MyPE = 0
      end subroutine initParent
!=====================================================================
!=====================================================================
!=====================================================================
      function isParent() result(p)
         implicit none
         logical :: p
         p = Parent
      end function isParent
!=====================================================================
!=====================================================================
!=====================================================================
         function isChild() result(p)
            implicit none
            logical :: p
            p = Child
         end function isChild
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine startChildTasks(num_hosts,host_name,exec_name,    &
     &                              num_tasks)
            implicit none
            integer, intent(in) :: num_hosts
            character (len=*), intent(in) :: host_name(:)
            character (len=*), intent(in) :: exec_name(:)
            integer, intent(in) :: num_tasks(:)
!
            integer :: i, n
            integer, allocatable :: nhp(:)
!
            if (.not.Parent) then
               call ErrorHandler('startChildTasks',                   &
     &                           'calling process is not a parent')
            else if (ChildTasksStarted) then
               call ErrorHandler('startChildTasks',                   &
     &                 'startChildTasks has already been called once')
            endif
            NumHosts = num_hosts
            NumPEs = 0
            do i=1, num_hosts
               NumPEs = NumPEs + num_tasks(i)
            enddo
            allocate( TaskID(0:NumPEs) )
            TaskID(NumPEs) = ParentTaskID
            ParentPE = NumPEs
            MyPE = ParentPE
            if (NumPEs == 0) then
!              -------------------------------------------------------
               call WarningHandler('startChildTasks',                 &
     &                             'No child tasks started')
!              -------------------------------------------------------
               return
            endif
!
            allocate( Hosts(num_hosts) )
!
            n = 0
            do i=1,num_hosts
               Hosts(i)%NumPEs = num_tasks(i)
               Hosts(i)%executable = exec_name(i)
               Hosts(i)%address = host_name(i)
               allocate( Hosts(i)%PEs(num_tasks(i)) )
            enddo
!
            ChildTasksStarted = .true.
!
         end subroutine startChildTasks
!=====================================================================
!=====================================================================
!=====================================================================
         function getHostPEs(i) result(p)
            implicit none
            integer, intent(in) :: i
            integer, pointer :: p(:)
            integer :: n
!
            if (.not.Parent) then
               call ErrorHandler('getHostPEs',                        &
     &                           'calling process is not a parent')
            else if (i < 1 .or. i > NumHosts) then
               call ErrorHandler('getHostPEs','invalid host index',i)
            endif
            n = Hosts(i)%NumPEs
            p => Hosts(i)%PEs(1:n)
         end function getHostPEs
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine initMPP_default_comm()
            implicit none
#ifdef MPI
            call initMPP_comm(MPI_COMM_WORLD)
#else
            call initMPP_comm(0)
#endif
         end subroutine initMPP_default_comm

         subroutine initMPP_comm(comm)
            implicit none
            integer (kind=IntKind) :: comm
!
            communicator = comm
#ifdef MPI
            call MPI_init(info)
#endif
!           ----------------------------------------------------------
            call calNumPEs(NumPEs)
            MyPE = getMyPEIndex()
!           ----------------------------------------------------------
            received = .false.
            packing = .false.
            sourceNode=0
            MaxWaits = MaxWaits_param
            allocate(MsgID(MaxWaits), MsgBuffer(MaxWaits),            &
                     MsgSRC(MaxWaits), WaitType(MaxWaits))
            MsgID(1:MaxWaits)=0
            MsgBuffer(1:MaxWaits)=0
            MsgSRC(1:MaxWaits)=-1
            WaitType(1:MaxWaits)='NONE'
            NumWaits=0; NumSends=0; NumRecvs=0; NumPacks=0
            bpos=0; bsize=0
#ifdef MPI
            call mpi_type_size(MPI_INTEGER,INTEGER_BYTES,info)
            call mpi_type_size(MPI_CHARACTER,CHARACTER_BYTES,info)
            call mpi_type_size(MPI_DOUBLE_PRECISION,REAL_BYTES,info)
            call mpi_type_size(MPI_DOUBLE_COMPLEX,COMPLEX_BYTES,info)
!
            if (MyPE==0) then
               print *,' '
               print *,'            ', &
     &      '********************************************************'
               print *,'            ', &
     &      '*                                                      *'
               print *,'            ', &
     &      '*              Fortran 90 Interface Module             *'
               print *,'            ', &
     &      '*                                                      *'
               print *,'            ', &
     &      '*                         For                          *'
               print *,'            ', &
     &      '*                                                      *'
               print *,'            ', &
     &      '*              MPI Message Passing Library             *'
               print *,'            ', &
     &      '*                                                      *'
               print *,'            ', &
     &      '*                   Version No. 2.0k                   *'
               print *,'            ', &
     &      '*                                                      *'
               print *,'            ', &
     &      '********************************************************'
               print *,' '
               print *,'CHARACTER_BYTES = ',CHARACTER_BYTES
               print *,'INTEGER_BYTES   = ',INTEGER_BYTES  
               print *,'REAL_BYTES      = ',REAL_BYTES      
               print *,'COMPLEX_BYTES   = ',COMPLEX_BYTES      
               print *,'ADDRESS_BYTES   = ',MPI_ADDRESS_KIND      
               print *,' '
               print *,'Number of Procs = ',NumPEs
               print *,' '
            endif
#endif
            MyPE_comm = MyPE
            NumPEs_comm = NumPEs
         end subroutine initMPP_comm
!=====================================================================
!=====================================================================
!=====================================================================
         function getParent() result(n)
            implicit none
            integer :: n
            n = ParentPE
         end function getParent
!=====================================================================
!=====================================================================
!=====================================================================
         function getNumPEs() result(n)
            implicit none
            integer :: n
            n = NumPEs_comm
         end function getNumPEs
!=====================================================================
!=====================================================================
!=====================================================================
         function getMyPE() result(n)
            implicit none
            integer :: n
            n = MyPE_comm
         end function getMyPE
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine killParent()
            implicit none
!
            if (.not.Child .and. .not.Parent) then
               call ErrorHandler('killParent','not a Child or Parent')
            else if (ParentKilled) then
               call WarningHandler('KillParent',                      &
     &                             'Parent cannot be killed twice')   
               return
            endif
            ParentKilled = .true.
         end subroutine killParent
!=====================================================================
!=====================================================================
!=====================================================================
         function isParentAlive() result(p)
            implicit none
            logical :: p
            if (Child .and. .not.ParentKilled) then
               p = .true.
            else if (Parent .and. .not.ParentKilled) then
               p = .true.
            else
               p = .false.
            endif
         end function isParentAlive
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine endParent()
            implicit none
            if (.not.ParentKilled) then
               call ErrorHandler('endParent',                         &
     &                           'Need to call killParent() first')
            endif
            print *,'Parent is leaving!'
         end subroutine endParent
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine endMPP()
            call syncAllPEs()
#ifdef MPI
!           ----------------------------------------------------------
            call MPI_finalize(info)
!           ----------------------------------------------------------
#endif
            deallocate(MsgID, MsgBuffer, MsgSRC, WaitType)
         end subroutine endMPP
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine calNumPEs(n)
            implicit none
            integer (kind=IntKind), intent(out) :: n
            integer :: MPP
!
            MPP = 1
!
#ifdef MPI
!           ----------------------------------------------------------
            call MPI_comm_size(communicator, MPP, info)
!           ----------------------------------------------------------
#endif
            n = MPP
         end subroutine calNumPEs
!=====================================================================
!=====================================================================
!=====================================================================
         function getMyPEIndex() result(n)
            implicit none
            integer (kind=IntKind) :: n
!
#ifdef MPI
!           ----------------------------------------------------------
            call MPI_comm_rank(communicator, n, info)
!           ----------------------------------------------------------
#else
            n=0
#endif
         end function getMyPEIndex
!=====================================================================
!=====================================================================
!=====================================================================
         function getSource() result(n)
            implicit none
            integer (kind=IntKind) :: n
            n = sourceNode
         end function getSource
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine syncAllPEs()
            if (Parent) then
               call WarningHandler('syncAllPEs',                      &
     &                             'Parent cannot sync with children')
               return
            endif
#ifdef MPI
            call MPI_barrier(communicator,info)
#endif
         end subroutine syncAllPEs
!
!=====================================================================
!=====================================================================
!!!!!    interface sendMessage
!!!!!       module procedure send_int0, send_real0, send_cmplx0
!!!!!       module procedure send_int1, send_real1, send_cmplx1
!!!!!       module procedure send_int2, send_real2, send_cmplx2
!!!!!       module procedure send_str0, send_str1,  send_str2,
!!!!!    end interface
!=====================================================================
!=====================================================================
         subroutine send_int0(a,msgtyp,dest)
            implicit none
            integer (kind=IntKind), intent(in) :: a
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            if (dest == MyPE_comm) then
!              ai=a
               return
            else if (dest.ge.NumPEs_comm) then
               print *,'ERROR at MyPE: bad target', MyPE_comm, dest
               stop
            endif
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_send(a,1,MPI_INTEGER,i,msgtyp,          &
     &                             communicator,info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_send(a,1,MPI_INTEGER,dest,msgtyp,             &
     &                       communicator,info)
!              -------------------------------------------------------
            endif
#endif
         end subroutine send_int0
!
!=====================================================================
         subroutine send_int1(a,n,msgtyp,dest)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            integer (kind=IntKind), intent(in) :: a(n)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            if (dest == MyPE_comm) then
!              allocate( ai1(n) )
!              ai1(1:n)=a(1:n)
               return
            else if (dest.ge.NumPEs_comm) then
               print *,'ERROR at MyPE: bad target', MyPE_comm, dest
               stop
            endif
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_send(a,n,MPI_INTEGER,i,msgtyp,          &
     &                             communicator,info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_send(a,n,MPI_INTEGER,dest,msgtyp,             &
     &                             communicator,info)
!              -------------------------------------------------------
            endif
#endif
         end subroutine send_int1
!
!=====================================================================
         subroutine send_int2(a,n,m,msgtyp,dest)
            implicit none
            integer (kind=IntKind), intent(in) :: n, m
            integer (kind=IntKind), intent(in) :: a(n,m)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            if (dest == MyPE_comm) then
!              allocate( ai2(n,m) )
!              ai2(1:n,1:m)=a(1:n,1:m)
               return
            else if (dest.ge.NumPEs_comm) then
               print *,'ERROR at MyPE: bad target', MyPE_comm, dest
               stop
            endif
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_send(a,n*m,MPI_INTEGER,i,msgtyp,        &
     &                             communicator,info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_send(a,n*m,MPI_INTEGER,dest,msgtyp,        &
     &                             communicator,info)
!              -------------------------------------------------------
            endif
#endif
         end subroutine send_int2
!
!=====================================================================
         subroutine send_real0(a,msgtyp,dest)
            implicit none
            real    (kind=RealKind), intent(in) :: a
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            if (dest == MyPE_comm) then
!              ar=a
               return
            else if (dest.ge.NumPEs_comm) then
               print *,'ERROR at MyPE: bad target', MyPE_comm, dest
               stop
            endif
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_send(a,1,MPI_DOUBLE_PRECISION,i,msgtyp, &
     &                             communicator,info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_send(a,1,MPI_DOUBLE_PRECISION,dest,msgtyp,    &
     &                             communicator,info)
!              -------------------------------------------------------
            endif
#endif
         end subroutine send_real0
!
!=====================================================================
         subroutine send_real1(a,n,msgtyp,dest)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            real    (kind=RealKind), intent(in) :: a(n)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            if (dest == MyPE_comm) then
!              allocate( ar1(n) )
!              ar1(1:n)=a(1:n)
               return
            else if (dest.ge.NumPEs_comm) then
               print *,'ERROR at MyPE: bad target', MyPE_comm, dest
               stop
            endif
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_send(a,n,MPI_DOUBLE_PRECISION,i,msgtyp, &
     &                             communicator,info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_send(a,n,MPI_DOUBLE_PRECISION,dest,msgtyp,    &
     &                             communicator,info)
!              -------------------------------------------------------
            endif
#endif
         end subroutine send_real1
!
!=====================================================================
         subroutine send_real2(a,n,m,msgtyp,dest)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            real    (kind=RealKind), intent(in) :: a(n,m)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            if (dest == MyPE_comm) then
!              allocate( ar2(n,m) )
!              ar2(1:n,1:m)=a(1:n,1:m)
               return
            else if (dest.ge.NumPEs_comm) then
               print *,'ERROR at MyPE: bad target', MyPE_comm, dest
               stop
            endif
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_send(a,n*m,MPI_DOUBLE_PRECISION,i,msgtyp,   &
     &                             communicator,info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_send(a,n*m,MPI_DOUBLE_PRECISION,dest,msgtyp,   &
     &                             communicator,info)
!              -------------------------------------------------------
            endif
#endif
         end subroutine send_real2
!
!=====================================================================
         subroutine send_cmplx0(a,msgtyp,dest)
            implicit none
            complex (kind=CmplxKind), intent(in) :: a
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            if (dest == MyPE_comm) then
!              ac=a
               return
            else if (dest.ge.NumPEs_comm) then
               print *,'ERROR at MyPE: bad target', MyPE_comm, dest
               stop
            endif
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_send(a,1,MPI_DOUBLE_COMPLEX,i,msgtyp,   &
     &                             communicator,info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_send(a,1,MPI_DOUBLE_COMPLEX,dest,msgtyp,   &
     &                       communicator,info)
!              -------------------------------------------------------
            endif
#endif
         end subroutine send_cmplx0
!
!=====================================================================
         subroutine send_cmplx1(a,n,msgtyp,dest)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            complex (kind=CmplxKind), intent(in) :: a(n)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            if (dest == MyPE_comm) then
!              allocate( ac1(n) )
!              ac1(1:n)=a(1:n)
               return
            else if (dest.ge.NumPEs_comm) then
               print *,'ERROR at MyPE: bad target', MyPE_comm, dest
               stop
            endif
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_send(a,n,MPI_DOUBLE_COMPLEX,i,msgtyp,   &
     &                             communicator,info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_send(a,n,MPI_DOUBLE_COMPLEX,dest,msgtyp,   &
     &                       communicator,info)
!              -------------------------------------------------------
            endif
#endif
         end subroutine send_cmplx1
!
!=====================================================================
         subroutine send_cmplx2(a,n,m,msgtyp,dest)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            complex (kind=CmplxKind), intent(in) :: a(n,m)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            if (dest == MyPE_comm) then
!              allocate( ac2(n,m) )
!              ac2(1:n,1:m)=a(1:n,1:m)
               return
            else if (dest.ge.NumPEs_comm) then
               print *,'ERROR at MyPE: bad target', MyPE_comm, dest
               stop
            endif
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_send(a,n*m,MPI_DOUBLE_COMPLEX,i,msgtyp, &
     &                             communicator,info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_send(a,n*m,MPI_DOUBLE_COMPLEX,dest,msgtyp, &
     &                       communicator,info)
!              -------------------------------------------------------
            endif
#endif
         end subroutine send_cmplx2
!
!=====================================================================
         subroutine send_str0(a,msgtyp,dest)
            implicit none
            character (len=*), intent(in) :: a
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            if (dest == MyPE_comm) then
!              slen=len(a)
!              allocate( as(slen) )
!              do i=1,slen
!                 as(i)=a(i:i)
!              enddo
               return
            else if (dest.ge.NumPEs_comm) then
               print *,'ERROR at MyPE: bad target', MyPE_comm, dest
               stop
            endif
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_send(a,len(a),MPI_CHARACTER,i,msgtyp, &
     &                             communicator,info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_send(a,len(a),MPI_CHARACTER,dest,msgtyp, &
     &                       communicator,info)
!              -------------------------------------------------------
            endif
#endif
         end subroutine send_str0
!
!=====================================================================
         subroutine send_str1(a,n,msgtyp,dest)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            character (len=*), intent(in), target :: a(n)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            if (dest == MyPE_comm) then
!              slen=len(a(n))
!              allocate( as1(slen,n) )
!              do j=1,n
!                 do i=1,slen
!                    as1(i,j)=a(j)(i:i)
!                 enddo
!              enddo
               return
            else if (dest.ge.NumPEs_comm) then
               print *,'ERROR at MyPE: bad target', MyPE_comm, dest
               stop
            endif
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_send(a,len(a(1))*n,MPI_CHARACTER,i,msgtyp, &
     &                             communicator,info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_send(a,len(a(1))*n,MPI_CHARACTER,dest,msgtyp, &
     &                       communicator,info)
!              -------------------------------------------------------
            endif
#endif
         end subroutine send_str1
!
!=====================================================================
         subroutine send_str2(a,n,m,msgtyp,dest)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            character (len=*), intent(in) :: a(n,m)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            if (dest == MyPE_comm) then
!              slen=len(a(n,m))
!              allocate( as2(slen,n,m) )
!              do k=1,n
!                 do j=1,n
!                    do i=1,slen
!                       as2(i,j,k)=a(j,k)(i:i)
!                    enddo
!                 enddo
!              enddo
               return
            else if (dest.ge.NumPEs_comm) then
               print *,'ERROR at MyPE: bad target', MyPE_comm, dest
               stop
            endif
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_send(a,len(a(1,1))*n*m,MPI_CHARACTER,i, &
     &                             msgtyp,communicator,info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_send(a,len(a(1,1))*n*m,MPI_CHARACTER,dest, &
     &                       msgtyp,communicator,info)
!              -------------------------------------------------------
            endif
#endif
         end subroutine send_str2
!
!=====================================================================
!=====================================================================
!!!!!    interface msendMessage
!!!!!       module procedure msend_int0, msend_real0, msend_cmplx0
!!!!!       module procedure msend_int1, msend_real1, msend_cmplx1
!!!!!       module procedure msend_int2, msend_real2, msend_cmplx2
!!!!!       module procedure msend_str0, msend_str1,  msend_str2,
!!!!!    end interface
!=====================================================================
!=====================================================================
         subroutine msend_int0(a,msgtyp,targets,count)
            implicit none
            integer (kind=IntKind), intent(in) :: a
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: count
            integer, intent(in) :: targets(count)
!
            if (.not.isParentAlive()) then
               if (NumPEs_comm == 1) then
                  return
               else if (count >= NumPEs_comm) then
                  call ErrorHandler('msendMessage',                   &
     &                          'target count >= NumPEs',count,NumPEs_comm)
               endif
               do i=1,count
                  if (targets(i) >= NumPEs_comm) then
                     call ErrorHandler('msendMessage',             &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            else if (count > NumPEs_comm) then
               call ErrorHandler('msendMessage',                      &
                                 'target count > NumPEs',count,NumPEs_comm)
               do i=1,count
                  if (targets(i) > NumPEs_comm) then
                     call ErrorHandler('msendMessage',                &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            endif
#ifdef MPI
            do i=1,count
               if (targets(i).ne.MyPE_comm) then
!                 ----------------------------------------------------
                  call MPI_send(a,1,MPI_INTEGER,targets(i),msgtyp, &
     &                          communicator,info)
!                 ----------------------------------------------------
               endif
            enddo
#endif
         end subroutine msend_int0
!
!=====================================================================
         subroutine msend_int1(a,n,msgtyp,targets,count)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            integer (kind=IntKind), intent(in) :: a(n)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: count
            integer, intent(in) :: targets(count)
!
            if (.not.isParentAlive()) then
               if (NumPEs_comm == 1) then
                  return
               else if (count >= NumPEs_comm) then
                  call ErrorHandler('msendMessage',                   &
     &                          'target count >= NumPEs',count,NumPEs_comm)
               endif
               do i=1,count
                  if (targets(i) >= NumPEs_comm) then
                     call ErrorHandler('msendMessage',             &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            else if (count > NumPEs_comm) then
               call ErrorHandler('msendMessage',                      &
                                 'target count > NumPEs',count,NumPEs_comm)
               do i=1,count
                  if (targets(i) > NumPEs_comm) then
                     call ErrorHandler('msendMessage',                &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            endif
#ifdef MPI
            do i=1,count
               if (targets(i).ne.MyPE_comm) then
!                 ----------------------------------------------------
                  call MPI_send(a,n,MPI_INTEGER,targets(i),msgtyp, &
     &                          communicator,info)
!                 ----------------------------------------------------
               endif
            enddo
#endif
         end subroutine msend_int1
!
!=====================================================================
         subroutine msend_int2(a,n,m,msgtyp,targets,count)
            implicit none
            integer (kind=IntKind), intent(in) :: n, m
            integer (kind=IntKind), intent(in) :: a(n,m)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: count
            integer, intent(in) :: targets(count)
!
            if (.not.isParentAlive()) then
               if (NumPEs_comm == 1) then
                  return
               else if (count >= NumPEs_comm) then
                  call ErrorHandler('msendMessage',                   &
     &                          'target count >= NumPEs',count,NumPEs_comm)
               endif
               do i=1,count
                  if (targets(i) >= NumPEs_comm) then
                     call ErrorHandler('msendMessage',             &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            else if (count > NumPEs_comm) then
               call ErrorHandler('msendMessage',                      &
                                 'target count > NumPEs',count,NumPEs_comm)
               do i=1,count
                  if (targets(i) > NumPEs_comm) then
                     call ErrorHandler('msendMessage',                &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            endif
#ifdef MPI
            do i=1,count
               if (targets(i).ne.MyPE_comm) then
!                 ----------------------------------------------------
                  call MPI_send(a,n*m,MPI_INTEGER,targets(i),msgtyp, &
     &                          communicator,info)
!                 ----------------------------------------------------
               endif
            enddo
#endif
         end subroutine msend_int2
!
!=====================================================================
         subroutine msend_real0(a,msgtyp,targets,count)
            implicit none
            real    (kind=RealKind), intent(in) :: a
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: count
            integer, intent(in) :: targets(count)
!
            if (.not.isParentAlive()) then
               if (NumPEs_comm == 1) then
                  return
               else if (count >= NumPEs_comm) then
                  call ErrorHandler('msendMessage',                   &
     &                          'target count >= NumPEs',count,NumPEs_comm)
               endif
               do i=1,count
                  if (targets(i) >= NumPEs_comm) then
                     call ErrorHandler('msendMessage',             &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            else if (count > NumPEs_comm) then
               call ErrorHandler('msendMessage',                      &
                                 'target count > NumPEs',count,NumPEs_comm)
               do i=1,count
                  if (targets(i) > NumPEs_comm) then
                     call ErrorHandler('msendMessage',                &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            endif
#ifdef MPI
            do i=1,count
               if (targets(i).ne.MyPE_comm) then
!                 ----------------------------------------------------
                  call MPI_send(a,1,MPI_DOUBLE_PRECISION,targets(i),msgtyp, &
     &                          communicator,info)
!                 ----------------------------------------------------
               endif
            enddo
#endif
         end subroutine msend_real0
!
!=====================================================================
         subroutine msend_real1(a,n,msgtyp,targets,count)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            real    (kind=RealKind), intent(in) :: a(n)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: count
            integer, intent(in) :: targets(count)
!
            if (.not.isParentAlive()) then
               if (NumPEs_comm == 1) then
                  return
               else if (count >= NumPEs_comm) then
                  call ErrorHandler('msendMessage',                   &
     &                          'target count >= NumPEs',count,NumPEs_comm)
               endif
               do i=1,count
                  if (targets(i) >= NumPEs_comm) then
                     call ErrorHandler('msendMessage',             &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            else if (count > NumPEs_comm) then
               call ErrorHandler('msendMessage',                      &
                                 'target count > NumPEs',count,NumPEs_comm)
               do i=1,count
                  if (targets(i) > NumPEs_comm) then
                     call ErrorHandler('msendMessage',                &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            endif
#ifdef MPI
            do i=1,count
               if (targets(i).ne.MyPE_comm) then
!                 ----------------------------------------------------
                  call MPI_send(a,n,MPI_DOUBLE_PRECISION,targets(i),msgtyp, &
     &                          communicator,info)
!                 ----------------------------------------------------
               endif
            enddo
#endif
         end subroutine msend_real1
!
!=====================================================================
         subroutine msend_real2(a,n,m,msgtyp,targets,count)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            real    (kind=RealKind), intent(in) :: a(n,m)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: count
            integer, intent(in) :: targets(count)
!
            if (.not.isParentAlive()) then
               if (NumPEs_comm == 1) then
                  return
               else if (count >= NumPEs_comm) then
                  call ErrorHandler('msendMessage',                   &
     &                          'target count >= NumPEs',count,NumPEs_comm)
               endif
               do i=1,count
                  if (targets(i) >= NumPEs_comm) then
                     call ErrorHandler('msendMessage',             &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            else if (count > NumPEs_comm) then
               call ErrorHandler('msendMessage',                      &
                                 'target count > NumPEs',count,NumPEs_comm)
               do i=1,count
                  if (targets(i) > NumPEs_comm) then
                     call ErrorHandler('msendMessage',                &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            endif
#ifdef MPI
            do i=1,count
               if (targets(i).ne.MyPE_comm) then
!                 ----------------------------------------------------
                  call MPI_send(a,n*m,MPI_DOUBLE_PRECISION,targets(i),msgtyp, &
     &                          communicator,info)
!                 ----------------------------------------------------
               endif
            enddo
#endif
         end subroutine msend_real2
!
!=====================================================================
         subroutine msend_cmplx0(a,msgtyp,targets,count)
            implicit none
            complex (kind=CmplxKind), intent(in) :: a
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: count
            integer, intent(in) :: targets(count)
!
            if (.not.isParentAlive()) then
               if (NumPEs_comm == 1) then
                  return
               else if (count >= NumPEs_comm) then
                  call ErrorHandler('msendMessage',                   &
     &                          'target count >= NumPEs',count,NumPEs_comm)
               endif
               do i=1,count
                  if (targets(i) >= NumPEs_comm) then
                     call ErrorHandler('msendMessage',             &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            else if (count > NumPEs_comm) then
               call ErrorHandler('msendMessage',                      &
                                 'target count > NumPEs',count,NumPEs_comm)
               do i=1,count
                  if (targets(i) > NumPEs_comm) then
                     call ErrorHandler('msendMessage',                &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            endif
#ifdef MPI
            do i=1,count
               if (targets(i).ne.MyPE_comm) then
!                 ----------------------------------------------------
                  call MPI_send(a,1,MPI_DOUBLE_COMPLEX,targets(i), &
     &                          msgtyp,communicator,info)
!                 ----------------------------------------------------
               endif
            enddo
#endif
         end subroutine msend_cmplx0
!
!=====================================================================
         subroutine msend_cmplx1(a,n,msgtyp,targets,count)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            complex (kind=CmplxKind), intent(in) :: a(n)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: count
            integer, intent(in) :: targets(count)
!
            if (.not.isParentAlive()) then
               if (NumPEs_comm == 1) then
                  return
               else if (count >= NumPEs_comm) then
                  call ErrorHandler('msendMessage',                   &
     &                          'target count >= NumPEs',count,NumPEs_comm)
               endif
               do i=1,count
                  if (targets(i) >= NumPEs_comm) then
                     call ErrorHandler('msendMessage',             &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            else if (count > NumPEs_comm) then
               call ErrorHandler('msendMessage',                      &
                                 'target count > NumPEs',count,NumPEs_comm)
               do i=1,count
                  if (targets(i) > NumPEs_comm) then
                     call ErrorHandler('msendMessage',                &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            endif
#ifdef MPI
            do i=1,count
               if (targets(i).ne.MyPE_comm) then
!                 ----------------------------------------------------
                  call MPI_send(a,n,MPI_DOUBLE_COMPLEX,targets(i), &
     &                          msgtyp,communicator,info)
!                 ----------------------------------------------------
               endif
            enddo
#endif
         end subroutine msend_cmplx1
!
!=====================================================================
         subroutine msend_cmplx2(a,n,m,msgtyp,targets,count)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            complex (kind=CmplxKind), intent(in) :: a(n,m)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: count
            integer, intent(in) :: targets(count)
!
            if (.not.isParentAlive()) then
               if (NumPEs_comm == 1) then
                  return
               else if (count >= NumPEs_comm) then
                  call ErrorHandler('msendMessage',                   &
     &                          'target count >= NumPEs',count,NumPEs_comm)
               endif
               do i=1,count
                  if (targets(i) >= NumPEs_comm) then
                     call ErrorHandler('msendMessage',             &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            else if (count > NumPEs_comm) then
               call ErrorHandler('msendMessage',                      &
                                 'target count > NumPEs',count,NumPEs_comm)
               do i=1,count
                  if (targets(i) > NumPEs_comm) then
                     call ErrorHandler('msendMessage',                &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            endif
#ifdef MPI
            do i=1,count
               if (targets(i).ne.MyPE_comm) then
!                 ----------------------------------------------------
                  call MPI_send(a,n*m,MPI_DOUBLE_COMPLEX,targets(i), &
     &                          msgtyp,communicator,info)
!                 ----------------------------------------------------
               endif
            enddo
#endif
         end subroutine msend_cmplx2
!
!=====================================================================
         subroutine msend_str0(a,msgtyp,targets,count)
            implicit none
            character (len=*), intent(in) :: a
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: count
            integer, intent(in) :: targets(count)
!
            if (.not.isParentAlive()) then
               if (NumPEs_comm == 1) then
                  return
               else if (count >= NumPEs_comm) then
                  call ErrorHandler('msendMessage',                   &
     &                          'target count >= NumPEs',count,NumPEs_comm)
               endif
               do i=1,count
                  if (targets(i) >= NumPEs_comm) then
                     call ErrorHandler('msendMessage',             &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            else if (count > NumPEs_comm) then
               call ErrorHandler('msendMessage',                      &
                                 'target count > NumPEs',count,NumPEs_comm)
               do i=1,count
                  if (targets(i) > NumPEs_comm) then
                     call ErrorHandler('msendMessage',                &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            endif
#ifdef MPI
            do i=1,count
               if (targets(i).ne.MyPE_comm) then
!                 ----------------------------------------------------
                  call MPI_send(a,len(a),MPI_CHARACTER,targets(i), &
     &                          msgtyp,communicator,info)
!                 ----------------------------------------------------
               endif
            enddo
#endif
         end subroutine msend_str0
!
!=====================================================================
         subroutine msend_str1(a,n,msgtyp,targets,count)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            character (len=*), intent(in), target :: a(n)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: count
            integer, intent(in) :: targets(count)
!
            if (.not.isParentAlive()) then
               if (NumPEs_comm == 1) then
                  return
               else if (count >= NumPEs_comm) then
                  call ErrorHandler('msendMessage',                   &
     &                          'target count >= NumPEs',count,NumPEs_comm)
               endif
               do i=1,count
                  if (targets(i) >= NumPEs_comm) then
                     call ErrorHandler('msendMessage',             &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            else if (count > NumPEs_comm) then
               call ErrorHandler('msendMessage',                      &
                                 'target count > NumPEs',count,NumPEs_comm)
               do i=1,count
                  if (targets(i) > NumPEs_comm) then
                     call ErrorHandler('msendMessage',                &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            endif
#ifdef MPI
            do i=1,count
               if (targets(i).ne.MyPE_comm) then
!                 ----------------------------------------------------
                  call MPI_send(a,len(a(n))*n,MPI_CHARACTER,targets(i), &
     &                          msgtyp,communicator,info)
!                 ----------------------------------------------------
               endif
            enddo
#endif
         end subroutine msend_str1
!
!=====================================================================
         subroutine msend_str2(a,n,m,msgtyp,targets,count)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            character (len=*), intent(in) :: a(n,m)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: count
            integer, intent(in) :: targets(count)
!
            if (.not.isParentAlive()) then
               if (NumPEs_comm == 1) then
                  return
               else if (count >= NumPEs_comm) then
                  call ErrorHandler('msendMessage',                   &
     &                          'target count >= NumPEs',count,NumPEs_comm)
               endif
               do i=1,count
                  if (targets(i) >= NumPEs_comm) then
                     call ErrorHandler('msendMessage',             &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            else if (count > NumPEs_comm) then
               call ErrorHandler('msendMessage',                      &
                                 'target count > NumPEs',count,NumPEs_comm)
               do i=1,count
                  if (targets(i) > NumPEs_comm) then
                     call ErrorHandler('msendMessage',                &
     &                                 'invalid target',targets(i))
                  endif
               enddo
            endif
#ifdef MPI
            do i=1,count
               if (targets(i).ne.MyPE_comm) then
!                 ----------------------------------------------------
                  call MPI_send(a,len(a(n,m))*n*m,MPI_CHARACTER, &
     &                          targets(i),msgtyp,communicator,info)
!                 ----------------------------------------------------
               endif
            enddo
#endif
         end subroutine msend_str2
!
!=====================================================================
!=====================================================================
!!!!!    interface recvMessage
!!!!!       module procedure recv_int0, recv_real0, recv_cmplx0
!!!!!       module procedure recv_int1, recv_real1, recv_cmplx1
!!!!!       module procedure recv_int2, recv_real2, recv_cmplx2
!!!!!       module procedure recv_str0, recv_str1,  recv_str2,
!!!!!       module procedure recvLattice, recvLmaxInfo
!!!!!    end interface
!=====================================================================
!=====================================================================
         subroutine recv_int0(a,msgtyp,source)
            implicit none
            integer (kind=IntKind), intent(inout) :: a
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            if (source == MyPE_comm) then
!              a=ai
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('recvMessage','source > NumPEs',     &
     &                           source,NumPEs_comm)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('recvMessage','source = NumPEs',     &
     &                           source,NumPEs_comm)
            endif
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_recv(a,1,MPI_INTEGER,source, &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_recv(a,1,MPI_INTEGER,MPI_ANY_SOURCE, &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=status(MPI_SOURCE)
            endif
#endif
         end subroutine recv_int0
!
!=====================================================================
         subroutine recv_int1(a,n,msgtyp,source)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            integer (kind=IntKind), intent(inout) :: a(n)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            if (source == MyPE_comm) then
!              a(1:n)=ai1(1:n)
!              deallocate( ai1 )
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('recvMessage','source > NumPEs',     &
     &                           source,NumPEs_comm)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('recvMessage','source = NumPEs',     &
     &                           source,NumPEs_comm)
            endif
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_recv(a,n,MPI_INTEGER,source, &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_recv(a,n,MPI_INTEGER,MPI_ANY_SOURCE, &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=status(MPI_SOURCE)
            endif
#endif
         end subroutine recv_int1
!
!=====================================================================
         subroutine recv_int2(a,n,m,msgtyp,source)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            integer (kind=IntKind), intent(inout) :: a(n,m)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            if (source == MyPE_comm) then
!              a(1:n,1:m)=ai2(1:n,1:m)
!              deallocate( ai2 )
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('recvMessage','source > NumPEs',     &
     &                           source,NumPEs_comm)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('recvMessage','source = NumPEs',     &
     &                           source,NumPEs_comm)
            endif
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_recv(a,n*m,MPI_INTEGER,source, &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_recv(a,n*m,MPI_INTEGER,MPI_ANY_SOURCE, &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=status(MPI_SOURCE)
            endif
#endif
         end subroutine recv_int2
!
!=====================================================================
         subroutine recv_real0(a,msgtyp,source)
            implicit none
            real    (kind=RealKind), intent(inout) :: a
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            if (source == MyPE_comm) then
!              a=ar
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('recvMessage','source > NumPEs',     &
     &                           source,NumPEs_comm)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('recvMessage','source = NumPEs',     &
     &                           source,NumPEs_comm)
            endif
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_recv(a,1,MPI_DOUBLE_PRECISION,source,         &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_recv(a,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=status(MPI_SOURCE)
            endif
#endif
         end subroutine recv_real0
!
!=====================================================================
         subroutine recv_real1(a,n,msgtyp,source)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            real    (kind=RealKind), intent(inout) :: a(n)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            if (source == MyPE_comm) then
!              a(1:n)=ar1(1:n)
!              deallocate( ar1 )
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('recvMessage','source > NumPEs',     &
     &                           source,NumPEs_comm)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('recvMessage','source = NumPEs',     &
     &                           source,NumPEs_comm)
            endif
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_recv(a,n,MPI_DOUBLE_PRECISION,source,         &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_recv(a,n,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=status(MPI_SOURCE)
            endif
#endif
         end subroutine recv_real1
!
!=====================================================================
         subroutine recv_real2(a,n,m,msgtyp,source)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            real    (kind=RealKind), intent(inout) :: a(n,m)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            if (source == MyPE_comm) then
!              a(1:n,1:m)=ar2(1:n,1:m)
!              deallocate( ar2 )
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('recvMessage','source > NumPEs',     &
     &                           source,NumPEs_comm)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('recvMessage','source = NumPEs',     &
     &                           source,NumPEs_comm)
            endif
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_recv(a,n*m,MPI_DOUBLE_PRECISION,source, &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_recv(a,n*m,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,    &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=status(MPI_SOURCE)
            endif
#endif
         end subroutine recv_real2
!
!=====================================================================
         subroutine recv_cmplx0(a,msgtyp,source)
            implicit none
            complex (kind=CmplxKind), intent(inout) :: a
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            if (source == MyPE_comm) then
!              a=ac
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('recvMessage','source > NumPEs',     &
     &                           source,NumPEs_comm)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('recvMessage','source = NumPEs',     &
     &                           source,NumPEs_comm)
            endif
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_recv(a,1,MPI_DOUBLE_COMPLEX,source, &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_recv(a,1,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE, &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=status(MPI_SOURCE)
            endif
#endif
         end subroutine recv_cmplx0
!
!=====================================================================
         subroutine recv_cmplx1(a,n,msgtyp,source)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            complex (kind=CmplxKind), intent(inout) :: a(n)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            if (source == MyPE_comm) then
!              a(1:n)=ac1(1:n)
!              deallocate( ac1 )
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('recvMessage','source > NumPEs',     &
     &                           source,NumPEs_comm)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('recvMessage','source = NumPEs',     &
     &                           source,NumPEs_comm)
            endif
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_recv(a,n,MPI_DOUBLE_COMPLEX,source, &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_recv(a,n,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE, &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=status(MPI_SOURCE)
            endif
#endif
         end subroutine recv_cmplx1
!
!=====================================================================
         subroutine recv_cmplx2(a,n,m,msgtyp,source)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            complex (kind=CmplxKind), intent(inout) :: a(n,m)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            if (source == MyPE_comm) then
!              a(1:n,1:m)=ac2(1:n,1:m)
!              deallocate( ac2 )
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('recvMessage','source > NumPEs',     &
     &                           source,NumPEs_comm)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('recvMessage','source = NumPEs',     &
     &                           source,NumPEs_comm)
            endif
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_recv(a,n*m,MPI_DOUBLE_COMPLEX,source, &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_recv(a,n*m,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE, &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=status(MPI_SOURCE)
            endif
#endif
         end subroutine recv_cmplx2
!
!=====================================================================
         subroutine recv_str0(a,msgtyp,source)
            implicit none
            character (len=*), intent(inout) :: a
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            if (source == MyPE_comm) then
!              l=min(slen,len(a))
!              do i=1,l
!                 a(i:i)=as(i)
!              enddo
!              deallocate( as )
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('recvMessage','source > NumPEs',     &
     &                           source,NumPEs_comm)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('recvMessage','source = NumPEs',     &
     &                           source,NumPEs_comm)
            endif
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_recv(a,len(a),MPI_CHARACTER,source, &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_recv(a,len(a),MPI_CHARACTER,MPI_ANY_SOURCE, &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=status(MPI_SOURCE)
            endif
#endif
         end subroutine recv_str0
!
!=====================================================================
         subroutine recv_str1(a,n,msgtyp,source)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            character (len=*), intent(inout) :: a(n)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            if (source == MyPE_comm) then
!              l=min(slen,len(a(n)))
!              do j=1,n
!                 do i=1,l
!                    a(j)(i:i)=as1(i,j)
!                 enddo
!              enddo
!              deallocate( as1 )
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('recvMessage','source > NumPEs',     &
     &                           source,NumPEs_comm)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('recvMessage','source = NumPEs',     &
     &                           source,NumPEs_comm)
            endif
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_recv(a,len(a(1))*n,MPI_CHARACTER,source, &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_recv(a,len(a(1))*n,MPI_CHARACTER,MPI_ANY_SOURCE, &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=status(MPI_SOURCE)
            endif
#endif
         end subroutine recv_str1
!
!=====================================================================
         subroutine recv_str2(a,n,m,msgtyp,source)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            character (len=*), intent(inout) :: a(n,m)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            if (source == MyPE_comm) then
!              l=min(slen,len(a(n,m)))
!              do k=1,m
!                 do j=1,n
!                    do i=1,l
!                       a(j,k)(i:i)=as2(i,j,k)
!                    enddo
!                 enddo
!              enddo
!              deallocate( as2 )
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('recvMessage','source > NumPEs',     &
     &                           source,NumPEs_comm)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('recvMessage','source = NumPEs',     &
     &                           source,NumPEs_comm)
            endif
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_recv(a,len(a(1,1))*n*m,MPI_CHARACTER,source, &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_recv(a,len(a(1,1))*n*m,MPI_CHARACTER, &
     &                       MPI_ANY_SOURCE, &
     &                       msgtyp,communicator,status,info)
!              -------------------------------------------------------
               sourceNode=status(MPI_SOURCE)
            endif
#endif
         end subroutine recv_str2
!
!=====================================================================
!=====================================================================
!!!!!    interface bcastMessage
!!!!!       module procedure bcast_int0, bcast_real0, bcast_cmplx0
!!!!!       module procedure bcast_int1, bcast_real1, bcast_cmplx1
!!!!!       module procedure bcast_int2, bcast_real2, bcast_cmplx2
!!!!!       module procedure bcast_str0, bcast_str1,  bcast_str2
!!!!!    end interface
!=====================================================================
!=====================================================================
         subroutine bcast_int0(a,root)
            implicit none
            integer (kind=IntKind), intent(inout) :: a
            integer, intent(in) :: root
            integer, parameter :: msgtyp = 10000001
            if (Parent) then
               call ErrorHandler('bcastMessage',                      &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            else if (root.ge.NumPEs_comm .or. root.lt.0) then
               call ErrorHandler('bcastMessage','invalid root',root)
            endif
            sourceNode=root
#ifdef MPI
!           ----------------------------------------------------------
            call MPI_bcast(a,1,MPI_INTEGER,root,communicator,info)
!           ----------------------------------------------------------
#endif
         end subroutine bcast_int0
!
!=====================================================================
         subroutine bcast_int1(a,n,root)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            integer (kind=IntKind), intent(inout) :: a(n)
            integer, intent(in) :: root
            integer, parameter :: msgtyp = 10000001
            if (Parent) then
               call ErrorHandler('bcastMessage',                      &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            else if (root.ge.NumPEs_comm .or. root.lt.0) then
               call ErrorHandler('bcastMessage','invalid root',root)
            endif
            sourceNode=root
#ifdef MPI
!           ----------------------------------------------------------
            call MPI_bcast(a,n,MPI_INTEGER,root,communicator,info)
!           ----------------------------------------------------------
#endif
         end subroutine bcast_int1
!
!=====================================================================
         subroutine bcast_int2(a,n,m,root)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            integer (kind=IntKind), intent(inout) :: a(n,m)
            integer, intent(in) :: root
            integer, parameter :: msgtyp = 10000001
            if (Parent) then
               call ErrorHandler('bcastMessage',                      &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            else if (root.ge.NumPEs_comm .or. root.lt.0) then
               call ErrorHandler('bcastMessage','invalid root',root)
            endif
            sourceNode=root
#ifdef MPI
!           ----------------------------------------------------------
            call MPI_bcast(a,n*m,MPI_INTEGER,root,communicator,info)
!           ----------------------------------------------------------
#endif
         end subroutine bcast_int2
!
!=====================================================================
         subroutine bcast_real0(a,root)
            implicit none
            real    (kind=RealKind), intent(inout) :: a
            integer, intent(in) :: root
            integer, parameter :: msgtyp = 10000001
            if (Parent) then
               call ErrorHandler('bcastMessage',                      &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            else if (root.ge.NumPEs_comm .or. root.lt.0) then
               call ErrorHandler('bcastMessage','invalid root',root)
            endif
            sourceNode=root
#ifdef MPI
!           ----------------------------------------------------------
            call MPI_bcast(a,1,MPI_DOUBLE_PRECISION,root,communicator,info)
!           ----------------------------------------------------------
#endif
         end subroutine bcast_real0
!
!=====================================================================
         subroutine bcast_real1(a,n,root)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            real    (kind=RealKind), intent(inout) :: a(n)
            integer, intent(in) :: root
            integer, parameter :: msgtyp = 10000001
            if (Parent) then
               call ErrorHandler('bcastMessage',                      &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            else if (root.ge.NumPEs_comm .or. root.lt.0) then
               call ErrorHandler('bcastMessage','invalid root',root)
            endif
            sourceNode=root
#ifdef MPI
!           ----------------------------------------------------------
            call MPI_bcast(a,n,MPI_DOUBLE_PRECISION,root,communicator,info)
!           ----------------------------------------------------------
#endif
         end subroutine bcast_real1
!
!=====================================================================
         subroutine bcast_real2(a,n,m,root)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            real    (kind=RealKind), intent(inout) :: a(n,m)
            integer, intent(in) :: root
            integer, parameter :: msgtyp = 10000001
            if (Parent) then
               call ErrorHandler('bcastMessage',                      &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            else if (root.ge.NumPEs_comm .or. root.lt.0) then
               call ErrorHandler('bcastMessage','invalid root',root)
            endif
            sourceNode=root
#ifdef MPI
!           ----------------------------------------------------------
            call MPI_bcast(a,n*m,MPI_DOUBLE_PRECISION,root,communicator,info)
!           ----------------------------------------------------------
#endif
         end subroutine bcast_real2
!
!=====================================================================
         subroutine bcast_cmplx0(a,root)
            implicit none
            complex (kind=CmplxKind), intent(inout) :: a
            integer, intent(in) :: root
            integer, parameter :: msgtyp = 10000001
            if (Parent) then
               call ErrorHandler('bcastMessage',                      &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            else if (root.ge.NumPEs_comm .or. root.lt.0) then
               call ErrorHandler('bcastMessage','invalid root',root)
            endif
            sourceNode=root
#ifdef MPI
!           ----------------------------------------------------------
            call MPI_bcast(a,1,MPI_DOUBLE_COMPLEX,root,communicator, &
     &                     info)
!           ----------------------------------------------------------
#endif
         end subroutine bcast_cmplx0
!
!=====================================================================
         subroutine bcast_cmplx1(a,n,root)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            complex (kind=CmplxKind), intent(inout) :: a(n)
            integer, intent(in) :: root
            integer, parameter :: msgtyp = 10000001
            if (Parent) then
               call ErrorHandler('bcastMessage',                      &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            else if (root.ge.NumPEs_comm .or. root.lt.0) then
               call ErrorHandler('bcastMessage','invalid root',root)
            endif
            sourceNode=root
#ifdef MPI
!           ----------------------------------------------------------
            call MPI_bcast(a,n,MPI_DOUBLE_COMPLEX,root,communicator, &
     &                     info)
!           ----------------------------------------------------------
#endif
         end subroutine bcast_cmplx1
!
!=====================================================================
         subroutine bcast_cmplx2(a,n,m,root)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            complex (kind=CmplxKind), intent(inout) :: a(n,m)
            integer, intent(in) :: root
            integer, parameter :: msgtyp = 10000001
            if (Parent) then
               call ErrorHandler('bcastMessage',                      &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            else if (root.ge.NumPEs_comm .or. root.lt.0) then
               call ErrorHandler('bcastMessage','invalid root',root)
            endif
            sourceNode=root
#ifdef MPI
!           ----------------------------------------------------------
            call MPI_bcast(a,n*m,MPI_DOUBLE_COMPLEX,root,communicator, &
     &                     info)
!           ----------------------------------------------------------
#endif
         end subroutine bcast_cmplx2
!
!=====================================================================
         subroutine bcast_str0(a,root)
            implicit none
            character (len=*), intent(inout) :: a
            integer, intent(in) :: root
            integer, parameter :: msgtyp = 10000001
            if (Parent) then
               call ErrorHandler('bcastMessage',                      &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            else if (root.ge.NumPEs_comm .or. root.lt.0) then
               call ErrorHandler('bcastMessage','invalid root',root)
            endif
            sourceNode=root
#ifdef MPI
!           ----------------------------------------------------------
            call MPI_bcast(a,len(a),MPI_CHARACTER,root,communicator, &
     &                     info)
!           ----------------------------------------------------------
#endif
         end subroutine bcast_str0
!
!=====================================================================
         subroutine bcast_str1(a,n,root)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            character (len=*), intent(inout) :: a(n)
            integer, intent(in) :: root
            integer, parameter :: msgtyp = 10000001
            if (Parent) then
               call ErrorHandler('bcastMessage',                      &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            else if (root.ge.NumPEs_comm .or. root.lt.0) then
               call ErrorHandler('bcastMessage','invalid root',root)
            endif
            sourceNode=root
#ifdef MPI
!           ----------------------------------------------------------
            call MPI_bcast(a,len(a(1))*n,MPI_CHARACTER,root, &
     &                     communicator,info)
!           ----------------------------------------------------------
#endif
         end subroutine bcast_str1
!
!=====================================================================
         subroutine bcast_str2(a,n,m,root)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            character (len=*), intent(inout) :: a(n,m)
            integer, intent(in) :: root
            integer, parameter :: msgtyp = 10000001
            if (Parent) then
               call ErrorHandler('bcastMessage',                      &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            else if (root.ge.NumPEs_comm .or. root.lt.0) then
               call ErrorHandler('bcastMessage','invalid root',root)
            endif
            sourceNode=root
#ifdef MPI
!           ----------------------------------------------------------
            call MPI_bcast(a,len(a(1,1))*n*m,MPI_CHARACTER,root, &
     &                     communicator,info)
!           ----------------------------------------------------------
#endif
         end subroutine bcast_str2
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine glb_max_int0(a)
            implicit none
            integer (kind=IntKind), intent(inout) :: a
            integer (kind=IntKind) :: b
            integer, parameter :: msgtyp = 10000011
!
            if (Parent) then
               call ErrorHandler('GlobalMax',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,1,MPI_INTEGER,MPI_MAX, &
     &                         communicator,info)
!           ----------------------------------------------------------
            a=b
#endif
         end subroutine glb_max_int0
!=====================================================================
         subroutine glb_max_int1(a,n)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            integer (kind=IntKind), intent(inout) :: a(n)
            integer (kind=IntKind), allocatable :: b(:)
            integer, parameter :: msgtyp = 10000011
!
            if (Parent) then
               call ErrorHandler('GlobalMax',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            allocate(b(1:n))
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,n,MPI_INTEGER,MPI_MAX, &
     &                         communicator,info)
!           ----------------------------------------------------------
            a(1:n)=b(1:n)
            deallocate(b)
#endif
         end subroutine glb_max_int1
!=====================================================================
         subroutine glb_max_int2(a,n1,n2)
            implicit none
            integer (kind=IntKind), intent(in) :: n1,n2
            integer (kind=IntKind), intent(inout) :: a(n1,n2)
            integer (kind=IntKind), allocatable :: b(:,:)
            integer, parameter :: msgtyp = 10000011
!
            if (Parent) then
               call ErrorHandler('GlobalMax',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            allocate(b(1:n1,1:n2))
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,n1*n2,MPI_INTEGER,MPI_MAX, &
     &                     communicator,info)
!           ----------------------------------------------------------
            a(1:n1,1:n2)=b(1:n1,1:n2)
            deallocate(b)
#endif
         end subroutine glb_max_int2
!=====================================================================
         subroutine glb_max_real0(a)
            implicit none
            integer, parameter :: msgtyp = 10000011
            real (kind=RealKind), intent(inout) :: a
            real (kind=RealKind) :: b
!
            if (Parent) then
               call ErrorHandler('GlobalMax',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
     &                         communicator,info)
!           ----------------------------------------------------------
            a=b
#endif
         end subroutine glb_max_real0
!=====================================================================
         subroutine glb_max_real1(a,n)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            integer, parameter :: msgtyp = 10000011
            real (kind=RealKind), intent(inout) :: a(n)
            real (kind=RealKind), allocatable :: b(:)
!
            if (Parent) then
               call ErrorHandler('GlobalMax',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            allocate(b(1:n))
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,n,MPI_DOUBLE_PRECISION,MPI_MAX, &
     &                         communicator,info)
!           ----------------------------------------------------------
            a(1:n)=b(1:n)
            deallocate(b)
#endif
         end subroutine glb_max_real1
!=====================================================================
         subroutine glb_max_real2(a,n1,n2)
            implicit none
            integer (kind=IntKind), intent(in) :: n1,n2
            integer, parameter :: msgtyp = 10000011
            real (kind=RealKind), intent(inout) :: a(n1,n2)
            real (kind=RealKind), allocatable :: b(:,:)
!
            if (Parent) then
               call ErrorHandler('GlobalMax',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            allocate(b(1:n1,1:n2))
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,n1*n2,MPI_DOUBLE_PRECISION,MPI_MAX, &
     &                         communicator,info)
!           ----------------------------------------------------------
            a(1:n1,1:n2)=b(1:n1,1:n2)
            deallocate(b)
#endif
         end subroutine glb_max_real2
!=====================================================================
         subroutine glb_max_cmplx0(a)
            implicit none
            integer, parameter :: msgtyp = 10000011
            complex (kind=CmplxKind), intent(inout) :: a
            complex (kind=CmplxKind) :: b
!
            if (Parent) then
               call ErrorHandler('GlobalMax',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,1,MPI_DOUBLE_COMPLEX,MPI_MAX, &
     &                         communicator,info)
!           ----------------------------------------------------------
            a=b
#endif
         end subroutine glb_max_cmplx0
!=====================================================================
         subroutine glb_max_cmplx1(a,n)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            integer, parameter :: msgtyp = 10000011
            complex (kind=CmplxKind), intent(inout) :: a(n)
            complex (kind=CmplxKind), allocatable :: b(:)
!
            if (Parent) then
               call ErrorHandler('GlobalMax',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            allocate(b(1:n))
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,n,MPI_DOUBLE_COMPLEX,MPI_MAX, &
     &                         communicator,info)
!           ----------------------------------------------------------
            a(1:n)=b(1:n)
            deallocate(b)
#endif
         end subroutine glb_max_cmplx1
!=====================================================================
         subroutine glb_max_cmplx2(a,n1,n2)
            implicit none
            integer (kind=IntKind), intent(in) :: n1,n2
            integer, parameter :: msgtyp = 10000011
            complex (kind=CmplxKind), intent(inout) :: a(n1,n2)
            complex (kind=CmplxKind), allocatable :: b(:,:)
!
            if (Parent) then
               call ErrorHandler('GlobalMax',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            allocate(b(1:n1,1:n2))
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,n1*n2,MPI_DOUBLE_COMPLEX,MPI_MAX, &
     &                         communicator,info)
!           ----------------------------------------------------------
            a(1:n1,1:n2)=b(1:n1,1:n2)
            deallocate(b)
#endif
         end subroutine glb_max_cmplx2
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine glb_min_int0(a)
            implicit none
            integer (kind=IntKind), intent(inout) :: a
            integer (kind=IntKind) :: b
            integer, parameter :: msgtyp = 10000011
!
            if (Parent) then
               call ErrorHandler('GlobalMin',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,1,MPI_INTEGER,MPI_MIN, &
     &                         communicator,info)
!           ----------------------------------------------------------
            a=b
#endif
         end subroutine glb_min_int0
!=====================================================================
         subroutine glb_min_int1(a,n)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            integer (kind=IntKind), intent(inout) :: a(n)
            integer (kind=IntKind), allocatable :: b(:)
            integer, parameter :: msgtyp = 10000011
!
            if (Parent) then
               call ErrorHandler('GlobalMin',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            allocate(b(1:n))
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,n,MPI_INTEGER,MPI_MIN,          &
     &                         communicator,info)
!           ----------------------------------------------------------
            a(1:n)=b(1:n)
            deallocate(b)
#endif
         end subroutine glb_min_int1
!=====================================================================
         subroutine glb_min_int2(a,n1,n2)
            implicit none
            integer (kind=IntKind), intent(in) :: n1,n2
            integer (kind=IntKind), intent(inout) :: a(n1,n2)
            integer (kind=IntKind), allocatable :: b(:,:)
            integer, parameter :: msgtyp = 10000011
!
            if (Parent) then
               call ErrorHandler('GlobalMin',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            allocate(b(1:n1,1:n2))
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,n1*n2,MPI_INTEGER,MPI_MIN,         &
     &                         communicator,info)
!           ----------------------------------------------------------
            a(1:n1,1:n2)=b(1:n1,1:n2)
            deallocate(b)
#endif
         end subroutine glb_min_int2
!=====================================================================
         subroutine glb_min_real0(a)
            implicit none
            integer, parameter :: msgtyp = 10000011
            real (kind=RealKind), intent(inout) :: a
            real (kind=RealKind) :: b
!
            if (Parent) then
               call ErrorHandler('GlobalMin',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,1,MPI_DOUBLE_PRECISION,MPI_MIN,    &
     &                         communicator,info)
!           ----------------------------------------------------------
            a=b
#endif
         end subroutine glb_min_real0
!=====================================================================
         subroutine glb_min_real1(a,n)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            integer, parameter :: msgtyp = 10000011
            real (kind=RealKind), intent(inout) :: a(n)
            real (kind=RealKind), allocatable :: b(:)
!
            if (Parent) then
               call ErrorHandler('GlobalMin',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            allocate(b(1:n))
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,n,MPI_DOUBLE_PRECISION,MPI_MIN,    &
     &                         communicator,info)
!           ----------------------------------------------------------
            a(1:n)=b(1:n)
            deallocate(b)
#endif
         end subroutine glb_min_real1
!=====================================================================
         subroutine glb_min_real2(a,n1,n2)
            implicit none
            integer (kind=IntKind), intent(in) :: n1,n2
            integer, parameter :: msgtyp = 10000011
            real (kind=RealKind), intent(inout) :: a(n1,n2)
            real (kind=RealKind), allocatable :: b(:,:)
!
            if (Parent) then
               call ErrorHandler('GlobalMin',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            allocate(b(1:n1,1:n2))
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,n1*n2,MPI_DOUBLE_PRECISION,MPI_MIN, &
     &                         communicator,info)
!           ----------------------------------------------------------
            a(1:n1,1:n2)=b(1:n1,1:n2)
            deallocate(b)
#endif
         end subroutine glb_min_real2
!=====================================================================
         subroutine glb_min_cmplx0(a)
            implicit none
            integer, parameter :: msgtyp = 10000011
            complex (kind=CmplxKind), intent(inout) :: a
            complex (kind=CmplxKind) :: b
!
            if (Parent) then
               call ErrorHandler('GlobalMin',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,1,MPI_DOUBLE_COMPLEX,MPI_MIN,       &
     &                         communicator,info)
!           ----------------------------------------------------------
            a=b
#endif
         end subroutine glb_min_cmplx0
!=====================================================================
         subroutine glb_min_cmplx1(a,n)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            integer, parameter :: msgtyp = 10000011
            complex (kind=CmplxKind), intent(inout) :: a(n)
            complex (kind=CmplxKind), allocatable :: b(:)
!
            if (Parent) then
               call ErrorHandler('GlobalMin',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            allocate(b(1:n))
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,n,MPI_DOUBLE_COMPLEX,MPI_MIN,      &
     &                         communicator,info)
!           ----------------------------------------------------------
            a(1:n)=b(1:n)
            deallocate(b)
#endif
         end subroutine glb_min_cmplx1
!=====================================================================
         subroutine glb_min_cmplx2(a,n1,n2)
            implicit none
            integer (kind=IntKind), intent(in) :: n1,n2
            integer, parameter :: msgtyp = 10000011
            complex (kind=CmplxKind), intent(inout) :: a(n1,n2)
            complex (kind=CmplxKind), allocatable :: b(:,:)
!
            if (Parent) then
               call ErrorHandler('GlobalMin',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            allocate(b(1:n1,1:n2))
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,n1*n2,MPI_DOUBLE_COMPLEX,MPI_MIN,  &
     &                         communicator,info)
!           ----------------------------------------------------------
            a(1:n1,1:n2)=b(1:n1,1:n2)
            deallocate(b)
#endif
         end subroutine glb_min_cmplx2
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine glb_sum_int0(a)
            implicit none
            integer (kind=IntKind), intent(inout) :: a
            integer (kind=IntKind) :: b
            integer, parameter :: msgtyp = 10000011
!
            if (Parent) then
               call ErrorHandler('GlobalSum',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            b = 0
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,1,MPI_INTEGER,MPI_SUM,             &
     &                         communicator,info)
!           ----------------------------------------------------------
            a=b
#endif
         end subroutine glb_sum_int0
!=====================================================================
         subroutine glb_sum_int1(a,n)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            integer (kind=IntKind), intent(inout) :: a(n)
            integer (kind=IntKind), allocatable :: b(:)
            integer, parameter :: msgtyp = 10000011
!
            if (Parent) then
               call ErrorHandler('GlobalSum',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            allocate(b(1:n))
            b = 0
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,n,MPI_INTEGER,MPI_SUM,            &
     &                         communicator,info)
!           ----------------------------------------------------------
            a(1:n)=b(1:n)
            deallocate(b)
#endif
         end subroutine glb_sum_int1
!=====================================================================
         subroutine glb_sum_int2(a,n1,n2)
            implicit none
            integer (kind=IntKind), intent(in) :: n1,n2
            integer (kind=IntKind), intent(inout) :: a(n1,n2)
            integer (kind=IntKind), allocatable :: b(:,:)
            integer, parameter :: msgtyp = 10000011
!
            if (Parent) then
               call ErrorHandler('GlobalSum',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            allocate(b(1:n1,1:n2))
            b = 0
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,n1*n2,MPI_INTEGER,MPI_SUM,         &
     &                         communicator,info)
!           ----------------------------------------------------------
            a(1:n1,1:n2)=b(1:n1,1:n2)
            deallocate(b)
#endif
         end subroutine glb_sum_int2
!=====================================================================
!=====================================================================
         subroutine glb_sum_int3(a,n1,n2,n3)
            implicit none
            integer (kind=IntKind), intent(in) :: n1,n2,n3
            integer (kind=IntKind), intent(inout) :: a(n1,n2,n3)
            integer (kind=IntKind), allocatable :: b(:,:,:)
            integer, parameter :: msgtyp = 10000011
!
            if (Parent) then
               call ErrorHandler('GlobalSum',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            allocate(b(1:n1,1:n2,1:n3))
            b = 0
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,n1*n2*n3,MPI_INTEGER,MPI_SUM,         &
     &                         communicator,info)
!           ----------------------------------------------------------
            a(1:n1,1:n2,1:n3)=b(1:n1,1:n2,1:n3)
            deallocate(b)
#endif
         end subroutine glb_sum_int3
!=====================================================================
         subroutine glb_sum_real0(a)
            implicit none
            integer, parameter :: msgtyp = 10000011
            real (kind=RealKind), intent(inout) :: a
            real (kind=RealKind) :: b
!
            if (Parent) then
               call ErrorHandler('GlobalSum',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            b = 0.0d0
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,1,MPI_DOUBLE_PRECISION,MPI_SUM,    &
     &                         communicator,info)
!           ----------------------------------------------------------
            a=b
#endif
         end subroutine glb_sum_real0
!=====================================================================
         subroutine glb_sum_real1(a,n)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            integer, parameter :: msgtyp = 10000011
            real (kind=RealKind), intent(inout) :: a(n)
            real (kind=RealKind), allocatable :: b(:)
!
            if (Parent) then
               call ErrorHandler('GlobalSum',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            allocate(b(1:n))
            b = 0.0d0
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,n,MPI_DOUBLE_PRECISION,MPI_SUM,    &
     &                         communicator,info)
!           ----------------------------------------------------------
            a(1:n)=b(1:n)
            deallocate(b)
#endif
         end subroutine glb_sum_real1
!=====================================================================
         subroutine glb_sum_real2(a,n1,n2)
            implicit none
            integer (kind=IntKind), intent(in) :: n1,n2
            integer, parameter :: msgtyp = 10000011
            real (kind=RealKind), intent(inout) :: a(n1,n2)
            real (kind=RealKind), allocatable :: b(:,:)
!
            if (Parent) then
               call ErrorHandler('GlobalSum',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            allocate(b(1:n1,1:n2))
            b = 0.0d0
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,n1*n2,MPI_DOUBLE_PRECISION,MPI_SUM, &
     &                         communicator,info)
!           ----------------------------------------------------------
            a(1:n1,1:n2)=b(1:n1,1:n2)
            deallocate(b)
#endif
         end subroutine glb_sum_real2
!=====================================================================

!=====================================================================
         subroutine glb_sum_real3(a,n1,n2,n3)
            implicit none
            integer (kind=IntKind), intent(in) :: n1,n2,n3
            integer, parameter :: msgtyp = 10000011
            real (kind=RealKind), intent(inout) :: a(n1,n2,n3)
            real (kind=RealKind), allocatable :: b(:,:,:)
!
            if (Parent) then
               call ErrorHandler('GlobalSum',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            allocate(b(1:n1,1:n2,1:n3))
            b = 0.0d0
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,n1*n2*n3,MPI_DOUBLE_PRECISION,MPI_SUM, &
     &                         communicator,info)
!           ----------------------------------------------------------
            a(1:n1,1:n2,1:n3)=b(1:n1,1:n2,1:n3)
            deallocate(b)
#endif
         end subroutine glb_sum_real3
!=====================================================================



         subroutine glb_sum_cmplx0(a)
            implicit none
            integer, parameter :: msgtyp = 10000011
            complex (kind=CmplxKind), intent(inout) :: a
            complex (kind=CmplxKind) :: b
!
            if (Parent) then
               call ErrorHandler('GlobalSum',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            b =  cmplx(0.d0,0.d0,kind=CmplxKind)
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,1,MPI_DOUBLE_COMPLEX,MPI_SUM,      &
     &                         communicator,info)
!           ----------------------------------------------------------
            a=b
#endif
         end subroutine glb_sum_cmplx0
!=====================================================================
         subroutine glb_sum_cmplx1(a,n)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            integer, parameter :: msgtyp = 10000011
            complex (kind=CmplxKind), intent(inout) :: a(n)
            complex (kind=CmplxKind), allocatable :: b(:)
!
            if (Parent) then
               call ErrorHandler('GlobalSum',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            allocate(b(1:n))
            b =  cmplx(0.d0,0.d0,kind=CmplxKind)
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,n,MPI_DOUBLE_COMPLEX,MPI_SUM,      &
     &                         communicator,info)
!           ----------------------------------------------------------
            a(1:n)=b(1:n)
            deallocate(b)
#endif
         end subroutine glb_sum_cmplx1
!=====================================================================
         subroutine glb_sum_cmplx2(a,n1,n2)
            implicit none
            integer (kind=IntKind), intent(in) :: n1,n2
            integer, parameter :: msgtyp = 10000011
            complex (kind=CmplxKind), intent(inout) :: a(n1,n2)
            complex (kind=CmplxKind), allocatable :: b(:,:)
!
            if (Parent) then
               call ErrorHandler('GlobalSum',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            allocate(b(1:n1,1:n2))
            b =  cmplx(0.d0,0.d0,kind=CmplxKind)
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,n1*n2,MPI_DOUBLE_COMPLEX,MPI_SUM,  &
     &                         communicator,info)
!           ----------------------------------------------------------
            a(1:n1,1:n2)=b(1:n1,1:n2)
            deallocate(b)
#endif
         end subroutine glb_sum_cmplx2
!=====================================================================
         subroutine glb_sum_cmplx3(a,n1,n2,n3)
            implicit none
            integer (kind=IntKind), intent(in) :: n1,n2,n3
            integer, parameter :: msgtyp = 10000011
            complex (kind=CmplxKind), intent(inout) :: a(n1,n2,n3)
            complex (kind=CmplxKind), allocatable :: b(:,:,:)
!
            if (Parent) then
               call ErrorHandler('GlobalSum',                         &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
#ifdef MPI
            allocate(b(1:n1,1:n2,1:n3))
            b =  cmplx(0.d0,0.d0,kind=CmplxKind)
!           ----------------------------------------------------------
            call MPI_ALLREDUCE(a,b,n1*n2*n3,MPI_DOUBLE_COMPLEX,MPI_SUM, &
     &                         communicator,info)
!           ----------------------------------------------------------
            a(1:n1,1:n2,1:n3)=b(1:n1,1:n2,1:n3)
            deallocate(b)
#endif
         end subroutine glb_sum_cmplx3
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine glb_collect_int0(a)
            implicit none
            integer (kind=IntKind), intent(inout) :: a(NumPEs_comm)
!
            if (Parent) then
               call ErrorHandler('GlobalCollect',                     &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
!
#ifdef MPI
!           ----------------------------------------------------------
!           call MPI_allgather(a(MyPE_comm+1),1,MPI_INTEGER,a,1,      &
            call MPI_allgather(MPI_IN_PLACE,1,MPI_INTEGER,a,1,        &
     &                         MPI_INTEGER,communicator,info)
!           ----------------------------------------------------------
#endif
         end subroutine glb_collect_int0
!=====================================================================
         subroutine glb_collect_int1(a,n)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            integer (kind=IntKind), intent(inout) :: a(n,NumPEs_comm)
!
            if (Parent) then
               call ErrorHandler('GlobalCollect',                     &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
!
#ifdef MPI
!           ----------------------------------------------------------
!           call MPI_allgather(a(1,MyPE_comm+1),n,MPI_INTEGER,a,n,    &
            call MPI_allgather(MPI_IN_PLACE,n,MPI_INTEGER,a,n,        &
     &                         MPI_INTEGER,communicator,info)
!           ----------------------------------------------------------
#endif
         end subroutine glb_collect_int1
!=====================================================================
         subroutine glb_collect_int2(a,n1,n2)
            implicit none
            integer (kind=IntKind), intent(in) :: n1,n2
            integer (kind=IntKind), intent(inout) :: a(n1,n2,NumPEs_comm)
!
            if (Parent) then
               call ErrorHandler('GlobalCollect',                     &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
!
#ifdef MPI
!           ----------------------------------------------------------
!           call MPI_allgather(a(1,1,MyPE_comm+1),n1*n2,MPI_INTEGER,a,&
            call MPI_allgather(MPI_IN_PLACE,n1*n2,MPI_INTEGER,a,      &
     &                         n1*n2,MPI_INTEGER,communicator,info)
!           ----------------------------------------------------------
#endif
         end subroutine glb_collect_int2
!=====================================================================
         subroutine glb_collect_real0(a)
            implicit none
            real (kind=RealKind), intent(inout) :: a(NumPEs_comm)
!
            if (Parent) then
               call ErrorHandler('GlobalCollect',                     &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
!
#ifdef MPI
!           ----------------------------------------------------------
!           call MPI_allgather(a(MyPE_comm+1),1,MPI_DOUBLE_PRECISION,a,1, &
            call MPI_allgather(MPI_IN_PLACE,1,MPI_DOUBLE_PRECISION,a,1,  &
     &                         MPI_DOUBLE_PRECISION,communicator,info)
!           ----------------------------------------------------------
#endif
         end subroutine glb_collect_real0
!=====================================================================
         subroutine glb_collect_real1(a,n)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            real (kind=RealKind), intent(inout) :: a(n,NumPEs_comm)
            real (kind=RealKind), allocatable :: b(:,:)
!
            if (Parent) then
               call ErrorHandler('GlobalCollect',                     &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
!
#ifdef MPI
         if (n <= MPP_BUFFER_MEM/8) then
!           ----------------------------------------------------------
!           call MPI_allgather(a(1,MyPE_comm+1),n,MPI_DOUBLE_PRECISION,a,n,  &
            call MPI_allgather(MPI_IN_PLACE,n,MPI_DOUBLE_PRECISION,a,n,  &
     &                         MPI_DOUBLE_PRECISION,communicator,info)
!           ----------------------------------------------------------
         else
            allocate(b(n,NumPEs_comm))
!           ----------------------------------------------------------
            call MPI_allgather(a(1,MyPE_comm+1),n,MPI_DOUBLE_PRECISION,b, &
     &                         n,MPI_DOUBLE_PRECISION,communicator,info)
!           ----------------------------------------------------------
            a = b
            deallocate(b)
         endif
#endif
         end subroutine glb_collect_real1
!=====================================================================
         subroutine glb_collect_real2(a,n1,n2)
            implicit none
            integer (kind=IntKind), intent(in) :: n1,n2
            real (kind=RealKind), intent(inout) :: a(n1,n2,NumPEs_comm)
            real (kind=RealKind), allocatable :: b(:,:,:)
!
            if (Parent) then
               call ErrorHandler('GlobalCollect',                     &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
!
#ifdef MPI
         if (n1*n2 <= MPP_BUFFER_MEM/8) then
!           ----------------------------------------------------------
!           call MPI_allgather(a(1,1,MyPE_comm+1),n1*n2,MPI_DOUBLE_PRECISION,a, &
            call MPI_allgather(MPI_IN_PLACE,n1*n2,MPI_DOUBLE_PRECISION,a, &
     &                         n1*n2,MPI_DOUBLE_PRECISION,communicator,info)
!           ----------------------------------------------------------
         else
            allocate(b(n1,n2,NumPEs_comm))
!           ----------------------------------------------------------
            call MPI_allgather(a(1,1,MyPE_comm+1),n1*n2,MPI_DOUBLE_PRECISION,b, &
     &                         n1*n2,MPI_DOUBLE_PRECISION,communicator,info)
!           ----------------------------------------------------------
            a = b
            deallocate(b)
         endif
#endif
         end subroutine glb_collect_real2
!=====================================================================
         subroutine glb_collect_cmplx0(a)
            implicit none
            complex (kind=CmplxKind), intent(inout) :: a(NumPEs_comm)
!
            if (Parent) then
               call ErrorHandler('GlobalCollect',                     &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
!
#ifdef MPI
!           ----------------------------------------------------------
!           call MPI_allgather(a(MyPE_comm+1),1,                      &
            call MPI_allgather(MPI_IN_PLACE,1,                        &
     &                         MPI_DOUBLE_COMPLEX,a,1,                &
     &                         MPI_DOUBLE_COMPLEX,communicator,info)
!           ----------------------------------------------------------
#endif
         end subroutine glb_collect_cmplx0
!=====================================================================
         subroutine glb_collect_cmplx1(a,n)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            complex (kind=CmplxKind), intent(inout) :: a(n,NumPEs_comm)
            complex (kind=CmplxKind), allocatable :: b(:,:)
!
            if (Parent) then
               call ErrorHandler('GlobalCollect',                     &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
!
#ifdef MPI
         if (n <= MPP_BUFFER_MEM/16) then
!           ----------------------------------------------------------
!           call MPI_allgather(a(1,MyPE_comm+1),n,                    &
            call MPI_allgather(MPI_IN_PLACE,n,                        &
     &                         MPI_DOUBLE_COMPLEX,a,n,                &
     &                         MPI_DOUBLE_COMPLEX,communicator,info)
!           ----------------------------------------------------------
         else
            allocate(b(n,NumPEs_comm))
!           ----------------------------------------------------------
            call MPI_allgather(a(1,MyPE_comm+1),n,                    &
     &                         MPI_DOUBLE_COMPLEX,b,n,                &
     &                         MPI_DOUBLE_COMPLEX,communicator,info)
!           ----------------------------------------------------------
            a = b
            deallocate(b)
         endif
#endif
         end subroutine glb_collect_cmplx1
!=====================================================================
         subroutine glb_collect_cmplx2(a,n1,n2)
            implicit none
            integer (kind=IntKind), intent(in) :: n1,n2
            complex (kind=CmplxKind), intent(inout) :: a(n1,n2,NumPEs_comm)
            complex (kind=CmplxKind), allocatable :: b(:,:,:)
!
            if (Parent) then
               call ErrorHandler('GlobalCollect',                     &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
!
#ifdef MPI
         if (n1*n2 <= MPP_BUFFER_MEM/16) then
!           ----------------------------------------------------------
!           call MPI_allgather(a(1,1,MyPE_comm+1),n1*n2,              &
            call MPI_allgather(MPI_IN_PLACE,n1*n2,                    &
     &                         MPI_DOUBLE_COMPLEX,a,n1*n2,            &
     &                         MPI_DOUBLE_COMPLEX,communicator,info)
!           ----------------------------------------------------------
         else
            allocate(b(n1,n2,NumPEs_comm))
!           ----------------------------------------------------------
            call MPI_allgather(a(1,1,MyPE_comm+1),n1*n2,              &
     &                         MPI_DOUBLE_COMPLEX,b,n1*n2,            &
     &                         MPI_DOUBLE_COMPLEX,communicator,info)
!           ----------------------------------------------------------
            a = b
            deallocate(b)
         endif
#endif
         end subroutine glb_collect_cmplx2
!=====================================================================
         subroutine glb_collect_str0(a)
            implicit none
            character (len=*), intent(inout) :: a(NumPEs_comm)
!
            if (Parent) then
               call ErrorHandler('GlobalCollect',                     &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
!
#ifdef MPI
!           ----------------------------------------------------------
!           call MPI_allgather(a(MyPE_comm+1),len(a(1)),              &
            call MPI_allgather(MPI_IN_PLACE,len(a(1)),                &
     &                         MPI_CHARACTER,a,len(a(1)),             &
     &                         MPI_CHARACTER,communicator,info)
!           ----------------------------------------------------------
#endif
         end subroutine glb_collect_str0
!=====================================================================
         subroutine glb_collect_str1(a,n)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            character (len=*), intent(inout) :: a(n,NumPEs_comm)
!
            if (Parent) then
               call ErrorHandler('GlobalCollect',                     &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
!
#ifdef MPI
!           ----------------------------------------------------------
!           call MPI_allgather(a(1,MyPE_comm+1),len(a(1,1))*n,        &
            call MPI_allgather(MPI_IN_PLACE,len(a(1,1))*n,            &
     &                         MPI_CHARACTER,a,len(a(1,1))*n,         &
     &                         MPI_CHARACTER,communicator,info)
!           ----------------------------------------------------------
#endif
         end subroutine glb_collect_str1
!=====================================================================
         subroutine glb_collect_str2(a,n1,n2)
            implicit none
            integer (kind=IntKind), intent(in) :: n1,n2
            character (len=*), intent(inout) :: a(n1,n2,NumPEs_comm)
!
            if (Parent) then
               call ErrorHandler('GlobalCollect',                     &
     &                   'Parent is not allowed to call bcastMessage')
            else if (NumPEs_comm == 1) then
               return
            endif
!
#ifdef MPI
!           call MPI_allgather(a(1,1,MyPE_comm+1),len(a(1,1,1))*n1*n2,     &
            call MPI_allgather(MPI_IN_PLACE,len(a(1,1,1))*n1*n2,      &
     &                         MPI_CHARACTER,a,len(a(1,1,1))*n1*n2,   &
     &                         MPI_CHARACTER,communicator,info)
!           ----------------------------------------------------------
#endif
         end subroutine glb_collect_str2
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine pack_int0(a)
            implicit none
            integer (kind=IntKind), intent(in) :: a
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.packing) then
               packing=.true.
               NumPacks=0
            endif
            NumPacks=NumPacks+1
#ifdef MPI
            PackedBuffer(NumPacks)%size=1
            PackedBuffer(NumPacks)%DataType=MPI_INTEGER
            allocate( PackedBuffer(NumPacks)%iBuf(1) )
            PackedBuffer(NumPacks)%iBuf(1)=a
#endif
         end subroutine pack_int0
!=====================================================================
         subroutine pack_int1(a,n)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            integer (kind=IntKind), intent(in) :: a(n)
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.packing) then
               packing=.true.
               NumPacks=0
            endif
            NumPacks=NumPacks+1
#ifdef MPI
            PackedBuffer(NumPacks)%size=n
            PackedBuffer(NumPacks)%DataType=MPI_INTEGER
            allocate( PackedBuffer(NumPacks)%iBuf(n) )
            PackedBuffer(NumPacks)%iBuf(1:n)=a(1:n)
#endif
         end subroutine pack_int1
!=====================================================================
         subroutine pack_int2(a,n,m)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            integer (kind=IntKind), intent(in) :: a(n,m)
            integer i, j
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.packing) then
               packing=.true.
               NumPacks=0
            endif
            NumPacks=NumPacks+1
#ifdef MPI
            PackedBuffer(NumPacks)%size=n*m
            PackedBuffer(NumPacks)%DataType=MPI_INTEGER
            allocate( PackedBuffer(NumPacks)%iBuf(n*m) )
            do i=1,m
               j=(i-1)*n
               PackedBuffer(NumPacks)%iBuf(j+1:j+n)=a(1:n,i)
            enddo
#endif
         end subroutine pack_int2
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine pack_longint0(a)
            implicit none
            integer (kind=LongIntKind), intent(in) :: a
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.packing) then
               packing=.true.
               NumPacks=0
            endif
            NumPacks=NumPacks+1
#ifdef MPI
            PackedBuffer(NumPacks)%size=2
            PackedBuffer(NumPacks)%DataType=MPI_INTEGER
            allocate( PackedBuffer(NumPacks)%iBuf(2) )
            PackedBuffer(NumPacks)%iBuf(1)=int(a/1000000000,kind=IntKind)
            PackedBuffer(NumPacks)%iBuf(2)=mod(a,1000000000_LongIntKind)
#endif
         end subroutine pack_longint0
!=====================================================================
         subroutine pack_longint1(a,n)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            integer (kind=LongIntKind), intent(in) :: a(n)
            integer (kind=IntKind) :: i
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.packing) then
               packing=.true.
               NumPacks=0
            endif
            NumPacks=NumPacks+1
#ifdef MPI
            PackedBuffer(NumPacks)%size=2*n
            PackedBuffer(NumPacks)%DataType=MPI_INTEGER
            allocate( PackedBuffer(NumPacks)%iBuf(2*n) )
            do i = 1, n
               PackedBuffer(NumPacks)%iBuf(2*i-1)=int(a(i)/1000000000,kind=IntKind)
               PackedBuffer(NumPacks)%iBuf(2*i)=mod(a(i),1000000000_LongIntKind)
            enddo
#endif
         end subroutine pack_longint1
!=====================================================================
         subroutine pack_longint2(a,n,m)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            integer (kind=LongIntKind), intent(in) :: a(n,m)
            integer i, j, k
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.packing) then
               packing=.true.
               NumPacks=0
            endif
            NumPacks=NumPacks+1
#ifdef MPI
            PackedBuffer(NumPacks)%size=2*n*m
            PackedBuffer(NumPacks)%DataType=MPI_INTEGER
            allocate( PackedBuffer(NumPacks)%iBuf(2*n*m) )
            do i=1,m
               j=2*(i-1)*n
               do k = 1, n
                  PackedBuffer(NumPacks)%iBuf(j+2*k-1)=int(a(k,i)/1000000000,kind=IntKind)
                  PackedBuffer(NumPacks)%iBuf(j+2*k)=mod(a(k,i),1000000000_LongIntKind)
               enddo
            enddo
#endif
         end subroutine pack_longint2
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine pack_real0(a)
            implicit none
            real (kind=RealKind), intent(in) :: a
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.packing) then
               packing=.true.
               NumPacks=0
            endif
            NumPacks=NumPacks+1
#ifdef MPI
            PackedBuffer(NumPacks)%size=1
            PackedBuffer(NumPacks)%DataType=MPI_DOUBLE_PRECISION
            allocate( PackedBuffer(NumPacks)%rBuf(1) )
            PackedBuffer(NumPacks)%rBuf(1)=a
#endif
         end subroutine pack_real0
!=====================================================================
         subroutine pack_real1(a,n)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            real (kind=RealKind), intent(in) :: a(n)
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.packing) then
               packing=.true.
               NumPacks=0
            endif
            NumPacks=NumPacks+1
#ifdef MPI
            PackedBuffer(NumPacks)%size=n
            PackedBuffer(NumPacks)%DataType=MPI_DOUBLE_PRECISION
            allocate( PackedBuffer(NumPacks)%rBuf(n) )
            PackedBuffer(NumPacks)%rBuf(1:n)=a(1:n)
#endif
         end subroutine pack_real1
!=====================================================================
         subroutine pack_real2(a,n,m)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            real (kind=RealKind), intent(in) :: a(n,m)
            integer i, j
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.packing) then
               packing=.true.
               NumPacks=0
            endif
            NumPacks=NumPacks+1
#ifdef MPI
            PackedBuffer(NumPacks)%size=n*m
            PackedBuffer(NumPacks)%DataType=MPI_DOUBLE_PRECISION
            allocate( PackedBuffer(NumPacks)%rBuf(n*m) )
            do i=1,m
               j=(i-1)*n
               PackedBuffer(NumPacks)%rBuf(j+1:j+n)=a(1:n,i)
            enddo
#endif
         end subroutine pack_real2
!=====================================================================
         subroutine pack_cmplx0(a)
            implicit none
            complex (kind=CmplxKind), intent(in) :: a
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.packing) then
               packing=.true.
               NumPacks=0
            endif
            NumPacks=NumPacks+1
#ifdef MPI
            PackedBuffer(NumPacks)%size=1
            PackedBuffer(NumPacks)%DataType=MPI_DOUBLE_COMPLEX
            allocate( PackedBuffer(NumPacks)%cBuf(1) )
            PackedBuffer(NumPacks)%cBuf(1)=a
#endif
         end subroutine pack_cmplx0
!=====================================================================
         subroutine pack_cmplx1(a,n)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            complex (kind=CmplxKind), intent(in) :: a(n)
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.packing) then
               packing=.true.
               NumPacks=0
            endif
            NumPacks=NumPacks+1
#ifdef MPI
            PackedBuffer(NumPacks)%size=n
            PackedBuffer(NumPacks)%DataType=MPI_DOUBLE_COMPLEX
            allocate( PackedBuffer(NumPacks)%cBuf(n) )
            PackedBuffer(NumPacks)%cBuf(1:n)=a(1:n)
#endif
         end subroutine pack_cmplx1
!=====================================================================
         subroutine pack_cmplx2(a,n,m)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            complex (kind=CmplxKind), intent(in) :: a(n,m)
            integer i, j
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.packing) then
               packing=.true.
               NumPacks=0
            endif
            NumPacks=NumPacks+1
#ifdef MPI
            PackedBuffer(NumPacks)%size=n*m
            PackedBuffer(NumPacks)%DataType=MPI_DOUBLE_COMPLEX
            allocate( PackedBuffer(NumPacks)%cBuf(n*m) )
            do i=1,m
               j=(i-1)*n
               PackedBuffer(NumPacks)%cBuf(j+1:j+n)=a(1:n,i)
            enddo
#endif
         end subroutine pack_cmplx2
!=====================================================================
         subroutine pack_str0(a)
            implicit none
            character (len=*), intent(in) :: a
            integer (kind=IntKind) :: size
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.packing) then
               packing=.true.
               NumPacks=0
            endif
            NumPacks=NumPacks+1
#ifdef MPI
            size=len(a)
            PackedBuffer(NumPacks)%size=size
            PackedBuffer(NumPacks)%DataType=MPI_CHARACTER
            allocate( PackedBuffer(NumPacks)%sBuf(size) )
            do i=1,size
               PackedBuffer(NumPacks)%sBuf(i)=a(i:i)
            enddo
#endif
         end subroutine pack_str0
!=====================================================================
         subroutine pack_str1(a,n)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            character (len=*), intent(in) :: a(n)
            integer (kind=IntKind) :: size
            integer (kind=IntKind) :: i, k, l, la
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.packing) then
               packing=.true.
               NumPacks=0
            endif
            NumPacks=NumPacks+1
#ifdef MPI
            la=len(a(n))
            size=n*la
            PackedBuffer(NumPacks)%size=size
            PackedBuffer(NumPacks)%DataType=MPI_CHARACTER
            allocate( PackedBuffer(NumPacks)%sBuf(size) )
!
            k=0
            do i=1,n
               do l=1,la
                  k=k+1
                  PackedBuffer(NumPacks)%sBuf(k)=a(i)(l:l)
               enddo
            enddo
!           call strcopy(a,n,PackedBuffer(NumPacks)%sBuf,size)
#endif
         end subroutine pack_str1
!=====================================================================
         subroutine pack_str2(a,n,m)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            character (len=*), intent(in) :: a(n,m)
            integer (kind=IntKind) :: size
            integer (kind=IntKind) :: i, j, k, l, la
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.packing) then
               packing=.true.
               NumPacks=0
            endif
            NumPacks=NumPacks+1
#ifdef MPI
            la=len(a(n,m))
            size=n*m*la
            PackedBuffer(NumPacks)%size=size
            PackedBuffer(NumPacks)%DataType=MPI_CHARACTER
            allocate( PackedBuffer(NumPacks)%sBuf(size) )
!
            k=0
            do j=1,m
               do i=1,n
                  do l=1,la
                     k=k+1
                     PackedBuffer(NumPacks)%sBuf(k)=a(i,j)(l:l)
                  enddo
               enddo
            enddo
!           call strcopy(a,n*m,PackedBuffer(NumPacks)%sBuf,size)
#endif
         end subroutine pack_str2
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine unpack_int0(a)
            implicit none
            integer (kind=IntKind), intent(inout) :: a
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.received) then
               print *,'ERROR:: need to call recvPackage first!'
               stop 'error'
            endif
#ifdef MPI
            call MPI_unpack(bspace,bsize,bpos,a,1,MPI_INTEGER,        &
     &                      communicator,info)
            if (bpos==bsize) then
               received=.false.
               deallocate( bspace )
            endif
#endif
         end subroutine unpack_int0
!=====================================================================
         subroutine unpack_int1(a,n)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            integer (kind=IntKind), intent(inout) :: a(n)
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.received) then
               print *,'ERROR:: need to call recvPackage first!'
               stop 'error'
            endif
#ifdef MPI
            call MPI_unpack(bspace,bsize,bpos,a,n,MPI_INTEGER,        &
     &                      communicator,info)
            if (bpos==bsize) then
               received=.false.
               deallocate( bspace )
            endif
#endif
         end subroutine unpack_int1
!=====================================================================
         subroutine unpack_int2(a,n,m)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            integer (kind=IntKind), intent(inout) :: a(n,m)
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.received) then
               print *,'ERROR:: need to call recvPackage first!'
               stop 'error'
            endif
#ifdef MPI
            call MPI_unpack(bspace,bsize,bpos,a,n*m,MPI_INTEGER,      &
     &                      communicator,info)
            if (bpos==bsize) then
               received=.false.
               deallocate( bspace )
            endif
#endif
         end subroutine unpack_int2
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine unpack_longint0(a)
            implicit none
            integer (kind=LongIntKind), intent(inout) :: a
            integer (kind=IntKind) :: at(2)
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.received) then
               print *,'ERROR:: need to call recvPackage first!'
               stop 'error'
            endif
#ifdef MPI
            call MPI_unpack(bspace,bsize,bpos,at,2,MPI_INTEGER,        &
     &                      communicator,info)
            if (bpos==bsize) then
               received=.false.
               deallocate( bspace )
            endif
            a = at(1)*1000000000+at(2)
#endif
         end subroutine unpack_longint0
!=====================================================================
         subroutine unpack_longint1(a,n)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            integer (kind=LongIntKind), intent(inout) :: a(n)
            integer (kind=IntKind) :: at(2*n), i
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.received) then
               print *,'ERROR:: need to call recvPackage first!'
               stop 'error'
            endif
#ifdef MPI
            call MPI_unpack(bspace,bsize,bpos,at,2*n,MPI_INTEGER,        &
     &                      communicator,info)
            if (bpos==bsize) then
               received=.false.
               deallocate( bspace )
            endif
            do i = 1, n
               a(i) = at(2*i-1)*1000000000 + at(2*i)
            enddo
#endif
         end subroutine unpack_longint1
!=====================================================================
         subroutine unpack_longint2(a,n,m)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            integer (kind=LongIntKind), intent(inout) :: a(n,m)
            integer (kind=IntKind) :: at(2*n,m), i, j
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.received) then
               print *,'ERROR:: need to call recvPackage first!'
               stop 'error'
            endif
#ifdef MPI
            call MPI_unpack(bspace,bsize,bpos,at,2*n*m,MPI_INTEGER,      &
     &                      communicator,info)
            if (bpos==bsize) then
               received=.false.
               deallocate( bspace )
            endif
            do j = 1, m
               do i = 1, n
                  a(i,j) = at(2*i-1,j)*1000000000 + at(2*i,j)
               enddo
            enddo
#endif
         end subroutine unpack_longint2
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine unpack_real0(a)
            implicit none
            real (kind=RealKind), intent(inout) :: a
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.received) then
               print *,'ERROR:: need to call recvPackage first!'
               stop 'error'
            endif
#ifdef MPI
            call MPI_unpack(bspace,bsize,bpos,a,1,MPI_DOUBLE_PRECISION, &
     &                      communicator,info)
            if (bpos==bsize) then
               received=.false.
               deallocate( bspace )
            endif
#endif
         end subroutine unpack_real0
!=====================================================================
         subroutine unpack_real1(a,n)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            real (kind=RealKind), intent(inout) :: a(n)
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.received) then
               print *,'ERROR:: need to call recvPackage first!'
               stop 'error'
            endif
#ifdef MPI
            call MPI_unpack(bspace,bsize,bpos,a,n,MPI_DOUBLE_PRECISION, &
     &                      communicator,info)
            if (bpos==bsize) then
               received=.false.
               deallocate( bspace )
            endif
#endif
         end subroutine unpack_real1
!=====================================================================
         subroutine unpack_real2(a,n,m)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            real (kind=RealKind), intent(inout) :: a(n,m)
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.received) then
               print *,'ERROR:: need to call recvPackage first!'
               stop 'error'
            endif
#ifdef MPI
            call MPI_unpack(bspace,bsize,bpos,a,n*m,MPI_DOUBLE_PRECISION, &
     &                      communicator,info)
            if (bpos==bsize) then
               received=.false.
               deallocate( bspace )
            endif
#endif
         end subroutine unpack_real2
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine unpack_cmplx0(a)
            implicit none
            complex (kind=CmplxKind), intent(inout) :: a
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.received) then
               print *,'ERROR:: need to call recvPackage first!'
               stop 'error'
            endif
#ifdef MPI
            call MPI_unpack(bspace,bsize,bpos,a,1,MPI_DOUBLE_COMPLEX, &
     &                      communicator,info)
            if (bpos==bsize) then
               received=.false.
               deallocate( bspace )
            endif
#endif
         end subroutine unpack_cmplx0
!=====================================================================
         subroutine unpack_cmplx1(a,n)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            complex (kind=CmplxKind), intent(inout) :: a(n)
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.received) then
               print *,'ERROR:: need to call recvPackage first!'
               stop 'error'
            endif
#ifdef MPI
            call MPI_unpack(bspace,bsize,bpos,a,n,MPI_DOUBLE_COMPLEX, &
     &                      communicator,info)
            if (bpos==bsize) then
               received=.false.
               deallocate( bspace )
            endif
#endif
         end subroutine unpack_cmplx1
!=====================================================================
         subroutine unpack_cmplx2(a,n,m)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            complex (kind=CmplxKind), intent(inout) :: a(n,m)
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.received) then
               print *,'ERROR:: need to call recvPackage first!'
               stop 'error'
            endif
#ifdef MPI
            call MPI_unpack(bspace,bsize,bpos,a,n*m,MPI_DOUBLE_COMPLEX, &
     &                      communicator,info)
            if (bpos==bsize) then
               received=.false.
               deallocate( bspace )
            endif
#endif
         end subroutine unpack_cmplx2
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine unpack_str0(a)
            implicit none
            character (len=*), intent(inout) :: a
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.received) then
               print *,'ERROR:: need to call recvPackage first!'
               stop 'error'
            endif
#ifdef MPI
            call MPI_unpack(bspace,bsize,bpos,a,len(a),MPI_CHARACTER, &
     &                      communicator,info)
            if (bpos==bsize) then
               received=.false.
               deallocate( bspace )
            endif
#endif
         end subroutine unpack_str0
!=====================================================================
         subroutine unpack_str1(a,n)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            character (len=*), intent(inout) :: a(n)
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.received) then
               print *,'ERROR:: need to call recvPackage first!'
               stop 'error'
            endif
#ifdef MPI
            call MPI_unpack(bspace,bsize,bpos,a,len(a(n))*n,          &
     &                      MPI_CHARACTER,communicator,info)
            if (bpos==bsize) then
               received=.false.
               deallocate( bspace )
            endif
#endif
         end subroutine unpack_str1
!=====================================================================
         subroutine unpack_str2(a,n,m)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            character (len=*), intent(inout) :: a(n,m)
            if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (.not.received) then
               print *,'ERROR:: need to call recvPackage first!'
               stop 'error'
            endif
#ifdef MPI
            call MPI_unpack(bspace,bsize,bpos,a,len(a(n,m))*n*m,      &
     &                      MPI_CHARACTER,communicator,info)
            if (bpos==bsize) then
               received=.false.
               deallocate( bspace )
            endif
#endif
         end subroutine unpack_str2
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine sendPackage0(msgtyp,dest)
            implicit none
            integer (kind=IntKind), intent(in) :: msgtyp
            integer (kind=IntKind), intent(in) :: dest
#ifdef MPI
            integer i, m
#endif
            if (dest == MyPE_comm) then
               return
            else if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (dest > NumPEs_comm) then
               call ErrorHandler('sendPackage','invalid target',dest)
            else if (.not.isParentAlive() .and. dest == NumPEs_comm) then
               call ErrorHandler('sendPackage','invalid target',dest)
            else if (.not.packing) then
               call ErrorHandler('sendPackage',                       &
     &                           'need to call packMessage first!')
            endif
#ifdef MPI
            bsize=0
            do i=1,NumPacks
!              -------------------------------------------------------
               call MPI_pack_size(PackedBuffer(i)%size,               &
     &                            PackedBuffer(i)%DataType,           &
     &                            communicator,m,info)
!              -------------------------------------------------------
               bsize=bsize+m
            enddo
!           ----------------------------------------------------------
            allocate( bspace(bsize) )
!           ----------------------------------------------------------
            call packBuffer()
!           ----------------------------------------------------------
            if (dest==AllPEs) then
               do i=1,NumPEs_comm
                  if (i-1.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_send(bspace,bpos,MPI_PACKED,i-1,msgtyp, &
     &                             communicator,info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_send(bspace,bpos,MPI_PACKED,dest,msgtyp,      &
     &                       communicator,info)
!              -------------------------------------------------------
            endif
!           ----------------------------------------------------------
            deallocate( bspace )
            call zeroBuffer()
!           ----------------------------------------------------------
#endif
            packing=.false.
         end subroutine sendPackage0
!=====================================================================
         subroutine sendPackage1(msgtyp,dest,n)
            implicit none
            integer (kind=IntKind), intent(in) :: msgtyp
            integer (kind=IntKind), intent(in) :: n
            integer (kind=IntKind), intent(in) :: dest(:)
            integer i
#ifdef MPI
            integer m
#endif
            if (n < 1) then
               call ErrorHandler('sendPackage','invalid no. targets',n)
            else if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (n == 1 .and. dest(1) == MyPE_comm) then
               return
            else if (.not.packing) then
               call ErrorHandler('sendPackage',                       &
     &                           'need to call packMessage first!')
            endif
#ifdef MPI
            bsize=0
            do i=1,NumPacks
!              -------------------------------------------------------
               call MPI_pack_size(PackedBuffer(i)%size,               &
     &                            PackedBuffer(i)%DataType,           &
     &                            communicator,m,info)
!              -------------------------------------------------------
               bsize=bsize+m
            enddo
            allocate( bspace(bsize) )
!           ----------------------------------------------------------
            call packBuffer()
!           ----------------------------------------------------------
            do i=1,n
!              -------------------------------------------------------
               call MPI_send(bspace,bpos,MPI_PACKED,dest(i),msgtyp,   &
     &                       communicator,info)
!              -------------------------------------------------------
            enddo
!           ----------------------------------------------------------
            deallocate( bspace )
!           ----------------------------------------------------------
            call zeroBuffer()
!           ----------------------------------------------------------
#endif
            packing=.false.
         end subroutine sendPackage1
!=====================================================================
!=====================================================================
!=====================================================================
         function nbsendPackage(msgtyp,dest) result(id)
            implicit none
            integer (kind=IntKind), intent(in) :: msgtyp
            integer (kind=IntKind), intent(in) :: dest
            integer (kind=IntKind) :: id
#ifdef MPI
            integer i, m
#endif
            if (dest == MyPE_comm) then
               id=1
               return
            else if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (dest > NumPEs_comm) then
               call ErrorHandler('nbsendPackage','invalid target',dest)
            else if (.not.isParentAlive() .and. dest == NumPEs_comm) then
               call ErrorHandler('nbsendPackage','invalid target',dest)
            else if (.not.packing) then
               call ErrorHandler('nbsendPackage',                     &
     &                           'need to call packMessage first!')
            endif
!
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbsendPackage',                     &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumSends=NumSends+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='SDPG'
#ifdef MPI
            bsize=0
            do i=1,NumPacks
!              -------------------------------------------------------
               call MPI_pack_size(PackedBuffer(i)%size,               &
     &                            PackedBuffer(i)%DataType,           &
     &                            communicator,m,info)
!              -------------------------------------------------------
               bsize=bsize+m
            enddo
!           ----------------------------------------------------------
            allocate( bspace(bsize) )
!           ----------------------------------------------------------
            call packBuffer()
!           ----------------------------------------------------------
            if (dest==AllPEs) then
               do i=1,NumPEs_comm
                  if (i-1.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_isend(bspace,bpos,MPI_PACKED,i-1,msgtyp, &
     &                              communicator,MsgID(id),info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_isend(bspace,bpos,MPI_PACKED,dest,msgtyp,     &
     &                        communicator,MsgID(id),info)
!              -------------------------------------------------------
            endif
!           ----------------------------------------------------------
            deallocate( bspace )
            call zeroBuffer()
!           ----------------------------------------------------------
#endif
            packing=.false.
         end function nbsendPackage
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine recvPackage(msgtyp,source)
            implicit none
            integer (kind=IntKind), intent(in) :: msgtyp
            integer (kind=IntKind), intent(in) :: source
            if (source == MyPE_comm) then
               return
            else if (.not.isParentAlive() .and. NumPEs_comm.eq.1) then
               return
            else if (source.ge.NumPEs_comm) then
               print *,'ERROR:: invalid source node: ',source
               stop
            endif
#ifdef MPI
            if (received) then
               print *,'WARNING:: MyPE = ',MyPE_comm, &
     &            ', the previous message is not been totally unpacked'
            endif
!           ----------------------------------------------------------
            call MPI_probe(source,msgtyp,communicator,status,info)
!           ----------------------------------------------------------
            call MPI_get_count(status,MPI_PACKED,bsize,info)
!           ----------------------------------------------------------
            if( allocated(bspace) ) then
               deallocate(bspace)
            endif
!           ----------------------------------------------------------
            allocate( bspace(bsize) )
!           ----------------------------------------------------------
            call MPI_recv(bspace,bsize,MPI_PACKED,source,msgtyp,      &
     &                    communicator,status,info)
!           ----------------------------------------------------------
#endif
            bpos=0             ! Need to add bpos=0 here to initialize
                               ! bpos for unpack functions
            received=.true.
         end subroutine recvPackage
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine strcopy(a,m,b,n)
            integer, intent(in) :: m,n
            integer :: l,i,j
            character (len=*), intent(in) :: a(m)
            character (len=1), intent(out) :: b(n)
            l = len(a(1))
            do i=1,m
               j = (i-1)*l
               b(j+1:j+l)=a(i)(1:l)
            enddo
         end subroutine strcopy
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine packBuffer()
#ifdef MPI
            integer i
            bpos=0
            do i=1,NumPacks
               if (PackedBuffer(i)%DataType==MPI_DOUBLE_COMPLEX) then
!                 ----------------------------------------------------
                  call MPI_pack(PackedBuffer(i)%cBuf,                 &
     &                          PackedBuffer(i)%size,                 &
     &                          MPI_DOUBLE_COMPLEX,                   &
     &                          bspace,bsize,bpos,communicator,info)
!                 ----------------------------------------------------
               else if (PackedBuffer(i)%DataType==MPI_INTEGER) then
!                 ----------------------------------------------------
                  call MPI_pack(PackedBuffer(i)%iBuf,                 &
     &                          PackedBuffer(i)%size,                 &
     &                          MPI_INTEGER,                          &
     &                          bspace,bsize,bpos,communicator,info)
!                 ----------------------------------------------------
               else if (PackedBuffer(i)%DataType==MPI_DOUBLE_PRECISION) then
!                 ----------------------------------------------------
                  call MPI_pack(PackedBuffer(i)%rBuf,                 &
     &                          PackedBuffer(i)%size,                 &
     &                          MPI_DOUBLE_PRECISION,                 &
     &                          bspace,bsize,bpos,communicator,info)
!                 ----------------------------------------------------
               else if (PackedBuffer(i)%DataType==MPI_CHARACTER) then
!                 ----------------------------------------------------
                  call MPI_pack(PackedBuffer(i)%sBuf,                 &
     &                          PackedBuffer(i)%size,                 &
     &                          MPI_CHARACTER,                        &
     &                          bspace,bsize,bpos,communicator,info)
!                 ----------------------------------------------------
               else 
                  print *,'ERROR:: unknown data type in packBuffer!'
                  stop
               endif
            enddo
#endif
         end subroutine packBuffer
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine zeroBuffer()
#ifdef MPI
            do i=1,NumPacks
               if (PackedBuffer(i)%DataType==MPI_DOUBLE_COMPLEX) then
                  deallocate( PackedBuffer(i)%cBuf )
                  nullify( PackedBuffer(i)%cBuf )
               else if (PackedBuffer(i)%DataType==MPI_INTEGER) then
                  deallocate( PackedBuffer(i)%iBuf )
                  nullify( PackedBuffer(i)%iBuf )
               else if (PackedBuffer(i)%DataType==MPI_DOUBLE_PRECISION) then
                  deallocate( PackedBuffer(i)%rBuf )
                  nullify( PackedBuffer(i)%rBuf )
               else if (PackedBuffer(i)%DataType==MPI_CHARACTER) then
                  deallocate( PackedBuffer(i)%sBuf )
                  nullify( PackedBuffer(i)%sBuf )
               endif
               PackedBuffer(i)%size=0
               PackedBuffer(i)%DataType=-1
            enddo
#endif
         end subroutine zeroBuffer
!=====================================================================
!=====================================================================
!=====================================================================
         function nbsend_int0(a,msgtyp,dest) result(id)
            implicit none
            integer (kind=IntKind), intent(in) :: a
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            integer (kind=IntKind) :: id
            if (dest == MyPE_comm) then
!              ai=a
               id=1
               return
            else if (dest > NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            else if (.not.isParentAlive() .and. dest == NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbsend_int0',                      &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumSends=NumSends+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='SEND'
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_isend(a,1,MPI_INTEGER,i,msgtyp,         &
     &                              communicator,MsgID(id),info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_isend(a,1,MPI_INTEGER,dest,msgtyp,         &
     &                        communicator,MsgID(id),info)
!              -------------------------------------------------------
            endif
#endif
         end function nbsend_int0
!=====================================================================
!=====================================================================
!=====================================================================
         function nbsend_real0(a,msgtyp,dest) result(id)
            implicit none
            real    (kind=RealKind), intent(in) :: a
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            integer (kind=IntKind) :: id
            if (dest == MyPE_comm) then
!              ar=a
               id=1
               return
            else if (dest > NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            else if (.not.isParentAlive() .and. dest == NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbsend_real0',                      &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumSends=NumSends+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='SEND'
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_isend(a,1,MPI_DOUBLE_PRECISION,i,msgtyp,      &
     &                             communicator,MsgID(id),info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_isend(a,1,MPI_DOUBLE_PRECISION,dest,msgtyp,         &
     &                        communicator,MsgID(id),info)
!              -------------------------------------------------------
            endif
#endif
         end function nbsend_real0
!=====================================================================
!=====================================================================
!=====================================================================
         function nbsend_cmplx0(a,msgtyp,dest) result(id)
            implicit none
            complex (kind=CmplxKind), intent(in) :: a
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            integer (kind=IntKind) :: id
            if (dest == MyPE_comm) then
!              ac=a
               id=1
               return
            else if (dest > NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            else if (.not.isParentAlive() .and. dest == NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbsend_cmplx0',                     &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumSends=NumSends+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='SEND'
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_isend(a,1,MPI_DOUBLE_COMPLEX,i,msgtyp,  &
     &                              communicator,MsgID(id),info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_isend(a,1,MPI_DOUBLE_COMPLEX,dest,msgtyp,  &
     &                        communicator,MsgID(id),info)
!              -------------------------------------------------------
            endif
#endif
         end function nbsend_cmplx0
!=====================================================================
!=====================================================================
!=====================================================================
         function nbsend_int1(a,n,msgtyp,dest) result(id)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            integer (kind=IntKind), intent(in) :: a(n)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            integer (kind=IntKind) :: id
            if (dest == MyPE_comm) then
!              allocate( ai1(n) )
!              ai1(1:n)=a(1:n)
               id=1
               return
            else if (dest > NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            else if (.not.isParentAlive() .and. dest == NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbsend_int1',                     &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumSends=NumSends+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='SEND'
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_isend(a,n,MPI_INTEGER,i,msgtyp,         &
     &                              communicator,MsgID(id),info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_isend(a,n,MPI_INTEGER,dest,msgtyp,         &
     &                        communicator,MsgID(id),info)
!              -------------------------------------------------------
            endif
#endif
         end function nbsend_int1
!=====================================================================
!=====================================================================
!=====================================================================
         function nbsend_real1(a,n,msgtyp,dest) result(id)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            real    (kind=RealKind), intent(in) :: a(n)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            integer (kind=IntKind) :: id
            if (dest == MyPE_comm) then
!              allocate( ar1(n) )
!              ar1(1:n)=a(1:n)
               id=1
               return
            else if (dest > NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            else if (.not.isParentAlive() .and. dest == NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbsend_real1',                     &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumSends=NumSends+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='SEND'
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_isend(a,n,MPI_DOUBLE_PRECISION,i,msgtyp,     &
     &                              communicator,MsgID(id),info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_isend(a,n,MPI_DOUBLE_PRECISION,dest,msgtyp,        &
     &                        communicator,MsgID(id),info)
!              -------------------------------------------------------
            endif
#endif
         end function nbsend_real1
!=====================================================================
!=====================================================================
!=====================================================================
         function nbsend_cmplx1(a,n,msgtyp,dest) result(id)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            complex (kind=CmplxKind), intent(in) :: a(n)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            integer (kind=IntKind) :: id
            if (dest == MyPE_comm) then
!              allocate( ac1(n) )
!              ac1(1:n)=a(1:n)
               id=1
               return
            else if (dest > NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            else if (.not.isParentAlive() .and. dest == NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbsend_cmplx1',                    &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumSends=NumSends+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='SEND'
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_isend(a,n,MPI_DOUBLE_COMPLEX,i,msgtyp,  &
     &                              communicator,MsgID(id),info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_isend(a,n,MPI_DOUBLE_COMPLEX,dest,msgtyp,     &
     &                        communicator,MsgID(id),info)
!              -------------------------------------------------------
            endif
#endif
         end function nbsend_cmplx1
!=====================================================================
!=====================================================================
!=====================================================================
         function nbsend_int2(a,n,m,msgtyp,dest) result(id)
            implicit none
            integer (kind=IntKind), intent(in) :: n, m
            integer (kind=IntKind), intent(in) :: a(n,m)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            integer (kind=IntKind) :: id
            if (dest == MyPE_comm) then
!              allocate( ai2(n,m) )
!              ai2(1:n,1:m)=a(1:n,1:m)
               id=1
               return
            else if (dest > NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            else if (.not.isParentAlive() .and. dest == NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbsend_int2',                       &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumSends=NumSends+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='SEND'
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_isend(a,n*m,MPI_INTEGER,i,msgtyp,       &
     &                              communicator,MsgID(id),info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_isend(a,n*m,MPI_INTEGER,dest,msgtyp,       &
     &                        communicator,MsgID(id),info)
!              -------------------------------------------------------
            endif
#endif
         end function nbsend_int2
!=====================================================================
!=====================================================================
!=====================================================================
         function nbsend_real2(a,n,m,msgtyp,dest) result(id)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            real    (kind=RealKind), intent(in) :: a(n,m)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            integer (kind=IntKind) :: id
            if (dest == MyPE_comm) then
!              allocate( ar2(n,m) )
!              ar2(1:n,1:m)=a(1:n,1:m)
               id=1
               return
            else if (dest > NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            else if (.not.isParentAlive() .and. dest == NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbsend_real2',                      &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumSends=NumSends+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='SEND'
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_isend(a,n*m,MPI_DOUBLE_PRECISION,i,msgtyp,   &
     &                              communicator,MsgID(id),info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_isend(a,n*m,MPI_DOUBLE_PRECISION,dest,msgtyp,      &
     &                        communicator,MsgID(id),info)
!              -------------------------------------------------------
            endif
#endif
         end function nbsend_real2
!=====================================================================
!=====================================================================
!=====================================================================
         function nbsend_cmplx2(a,n,m,msgtyp,dest) result(id)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            complex (kind=CmplxKind), intent(in) :: a(n,m)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            integer (kind=IntKind) :: id
            if (dest == MyPE_comm) then
!              allocate( ac2(n,m) )
!              ac2(1:n,1:m)=a(1:n,1:m)
               id=1
               return
            else if (dest > NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            else if (.not.isParentAlive() .and. dest == NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbsend_cmplx2',                     &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumSends=NumSends+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='SEND'
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_isend(a,n*m,MPI_DOUBLE_COMPLEX,i,msgtyp, &
     &                              communicator,MsgID(id),info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_isend(a,n*m,MPI_DOUBLE_COMPLEX,dest,msgtyp, &
     &                        communicator,MsgID(id),info)
!              -------------------------------------------------------
            endif
#endif
         end function nbsend_cmplx2
!=====================================================================
!=====================================================================
!=====================================================================
         function nbsend_str0(a,msgtyp,dest) result(id)
            implicit none
            character (len=*), intent(in) :: a
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            integer (kind=IntKind) :: id
            if (dest == MyPE_comm) then
!              slen=len(a)
!              allocate( as(slen) )
!              do i=1,slen
!                 as(i)=a(i:i)
!              enddo
               id=1
               return
            else if (dest > NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            else if (.not.isParentAlive() .and. dest == NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbsend_str0',                       &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumSends=NumSends+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='SEND'
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_isend(a,len(a),MPI_CHARACTER,i,msgtyp,  &
     &                              communicator,MsgID(id),info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_isend(a,len(a),MPI_CHARACTER,dest,msgtyp,     &
     &                        communicator,MsgID(id),info)
!              -------------------------------------------------------
            endif
#endif
         end function nbsend_str0
!=====================================================================
!=====================================================================
!=====================================================================
         function nbsend_str1(a,n,msgtyp,dest) result(id)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            character (len=*), intent(in), target :: a(n)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            integer (kind=IntKind) :: id
            if (dest == MyPE_comm) then
!              slen=len(a(n))
!              allocate( as1(slen,n) )
!              do j=1,n
!                 do i=1,slen
!                    as1(i,j)=a(j)(i:i)
!                 enddo
!              enddo
               id=1
               return
            else if (dest > NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            else if (.not.isParentAlive() .and. dest == NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbsend_str1',                       &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumSends=NumSends+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='SEND'
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_isend(a,len(a(1))*n,MPI_CHARACTER,i,    &
     &                           msgtyp,communicator,MsgID(id),info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_isend(a,len(a(1))*n,MPI_CHARACTER,dest,msgtyp, &
     &                        communicator,MsgID(id),info)
!              -------------------------------------------------------
            endif
#endif
         end function nbsend_str1
!=====================================================================
!=====================================================================
!=====================================================================
         function nbsend_str2(a,n,m,msgtyp,dest) result(id)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            character (len=*), intent(in) :: a(n,m)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: dest
            integer (kind=IntKind) :: id
            if (dest == MyPE_comm) then
!              slen=len(a(n,m))
!              allocate( as2(slen,n,m) )
!              do k=1,n
!                 do j=1,n
!                    do i=1,slen
!                       as2(i,j,k)=a(j,k)(i:i)
!                    enddo
!                 enddo
!              enddo
               id=1
               return
            else if (dest > NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            else if (.not.isParentAlive() .and. dest == NumPEs_comm) then
               call ErrorHandler('nbsendMessage','invalid target',dest)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbsend_str2',                       &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumSends=NumSends+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='SEND'
#ifdef MPI
            if (dest == All) then
               do i=0,NumPEs_comm-1
                  if (i.ne.MyPE_comm) then
!                    -------------------------------------------------
                     call MPI_isend(a,len(a(1,1))*n*m,MPI_CHARACTER,i,&
     &                           msgtyp,communicator,MsgID(id),info)
!                    -------------------------------------------------
                  endif
               enddo
            else
!              -------------------------------------------------------
               call MPI_isend(a,len(a(1,1))*n*m,MPI_CHARACTER,dest,   &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
            endif
#endif
         end function nbsend_str2
!=====================================================================
!=====================================================================
!=====================================================================
         function nbrecv_int0(a,msgtyp,source) result(id)
            implicit none
            integer (kind=IntKind), intent(inout), target :: a
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            integer (kind=IntKind) :: id
            if (source == MyPE_comm) then
!              a=ai
               id=1
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbrecv_int0',                       &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumRecvs=NumRecvs+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='RECV'
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_irecv(a,1,MPI_INTEGER,source,    &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_irecv(a,1,MPI_INTEGER,MPI_ANY_SOURCE,    &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=-1
            endif
#endif
         end function nbrecv_int0
!=====================================================================
!=====================================================================
!=====================================================================
         function nbrecv_real0(a,msgtyp,source) result(id)
            implicit none
            real    (kind=RealKind), intent(inout), target :: a
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            integer (kind=IntKind) :: id
            if (source == MyPE_comm) then
!              a=ar
               id=1
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbrecv_real0',                      &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumRecvs=NumRecvs+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='RECV'
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_irecv(a,1,MPI_DOUBLE_PRECISION,source,    &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_irecv(a,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,    &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=-1
            endif
#endif
         end function nbrecv_real0
!=====================================================================
!=====================================================================
!=====================================================================
         function nbrecv_cmplx0(a,msgtyp,source) result(id)
            implicit none
            complex (kind=CmplxKind), intent(inout), target :: a
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            integer (kind=IntKind) :: id
            if (source == MyPE_comm) then
!              a=ac
               id=1
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbrecv_cmplx0',                     &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumRecvs=NumRecvs+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='RECV'
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_irecv(a,1,MPI_DOUBLE_COMPLEX,source,   &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_irecv(a,1,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE, &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=-1
            endif
#endif
         end function nbrecv_cmplx0
!=====================================================================
!=====================================================================
!=====================================================================
         function nbrecv_int1(a,n,msgtyp,source) result(id)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            integer (kind=IntKind), intent(inout), target :: a(n)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            integer (kind=IntKind) :: id
            if (source == MyPE_comm) then
!              a(1:n)=ai1(1:n)
!              deallocate( ai1 )
               id=1
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbrecv_int1',                       &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumRecvs=NumRecvs+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='RECV'
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_irecv(a,n,MPI_INTEGER,source, &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_irecv(a,n,MPI_INTEGER,MPI_ANY_SOURCE, &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=-1
            endif
#endif
         end function nbrecv_int1
!=====================================================================
!=====================================================================
!=====================================================================
         function nbrecv_real1(a,n,msgtyp,source) result(id)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            real    (kind=RealKind), intent(inout), target :: a(n)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            integer (kind=IntKind) :: id
            if (source == MyPE_comm) then
!              a(1:n)=ar1(1:n)
!              deallocate( ar1 )
               id=1
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbrecv_real1',                      &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumRecvs=NumRecvs+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='RECV'
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_irecv(a,n,MPI_DOUBLE_PRECISION,source, &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_irecv(a,n,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=-1
            endif
#endif
         end function nbrecv_real1
!=====================================================================
!=====================================================================
!=====================================================================
         function nbrecv_cmplx1(a,n,msgtyp,source) result(id)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            complex (kind=CmplxKind), intent(inout), target :: a(n)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            integer (kind=IntKind) :: id
            if (source == MyPE_comm) then
!              a(1:n)=ac1(1:n)
!              deallocate( ac1 )
               id=1
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbrecv_cmplx1',                     &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumRecvs=NumRecvs+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='RECV'
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_irecv(a,n,MPI_DOUBLE_COMPLEX,source, &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_irecv(a,n,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE, &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=-1
            endif
#endif
         end function nbrecv_cmplx1
!=====================================================================
!=====================================================================
!=====================================================================
         function nbrecv_int2(a,n,m,msgtyp,source) result(id)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            integer (kind=IntKind), intent(inout), target :: a(n,m)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            integer (kind=IntKind) :: id
            if (source == MyPE_comm) then
!              a(1:n,1:m)=ai2(1:n,1:m)
!              deallocate( ai2 )
               id=1
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbrecv_int2',                       &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumRecvs=NumRecvs+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='RECV'
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_irecv(a,n*m,MPI_INTEGER,source, &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_irecv(a,n*m,MPI_INTEGER,MPI_ANY_SOURCE, &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=-1
            endif
#endif
         end function nbrecv_int2
!=====================================================================
!=====================================================================
!=====================================================================
         function nbrecv_real2(a,n,m,msgtyp,source) result(id)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            real    (kind=RealKind), intent(inout), target :: a(n,m)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            integer (kind=IntKind) :: id
            if (source == MyPE_comm) then
!              a(1:n,1:m)=ar2(1:n,1:m)
!              deallocate( ar2 )
               id=1
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbrecv_real2',                      &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumRecvs=NumRecvs+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='RECV'
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_irecv(a,n*m,MPI_DOUBLE_PRECISION,source, &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_irecv(a,n*m,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,    &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=-1
            endif
#endif
         end function nbrecv_real2
!=====================================================================
!=====================================================================
!=====================================================================
         function nbrecv_cmplx2(a,n,m,msgtyp,source) result(id)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            complex (kind=CmplxKind), intent(inout), target :: a(n,m)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            integer (kind=IntKind) :: id
            if (source == MyPE_comm) then
!              a(1:n,1:m)=ac2(1:n,1:m)
!              deallocate( ac2 )
               id=1
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbrecv_cmplx2',                     &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumRecvs=NumRecvs+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='RECV'
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_irecv(a,n*m,MPI_DOUBLE_COMPLEX,source, &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_irecv(a,n*m,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE, &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=-1
            endif
#endif
         end function nbrecv_cmplx2
!=====================================================================
!=====================================================================
!=====================================================================
         function nbrecv_str0(a,msgtyp,source) result(id)
            implicit none
            character (len=*), intent(inout), target :: a
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            integer (kind=IntKind) :: id
            integer (kind=IntKind) :: la
            if (source == MyPE_comm) then
!              l=min(slen,len(a))
!              do i=1,l
!                 a(i:i)=as(i)
!              enddo
!              deallocate( as )
               id=1
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbrecv_str0',                      &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumRecvs=NumRecvs+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='RECV'
            la=len(a)
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_irecv(a,len(a),MPI_CHARACTER,source, &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_irecv(a,len(a),MPI_CHARACTER,MPI_ANY_SOURCE, &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=-1
            endif
#endif
         end function nbrecv_str0
!=====================================================================
!=====================================================================
!=====================================================================
         function nbrecv_str1(a,n,msgtyp,source) result(id)
            implicit none
            integer (kind=IntKind), intent(in) :: n
            character (len=*), intent(inout), target :: a(n)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            integer (kind=IntKind) :: id
            integer (kind=IntKind) :: la
            if (source == MyPE_comm) then
!              l=min(slen,len(a(n)))
!              do j=1,n
!                 do i=1,l
!                    a(j)(i:i)=as1(i,j)
!                 enddo
!              enddo
!              deallocate( as1 )
               id=1
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbrecv_str1',                      &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumRecvs=NumRecvs+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='RECV'
            la=len(a(n))
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_irecv(a,len(a(1))*n,MPI_CHARACTER,source, &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_irecv(a,len(a(1))*n,MPI_CHARACTER,              &
     &                        MPI_ANY_SOURCE, &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=-1
            endif
#endif
         end function nbrecv_str1
!=====================================================================
!=====================================================================
!=====================================================================
         function nbrecv_str2(a,n,m,msgtyp,source) result(id)
            implicit none
            integer (kind=IntKind), intent(in) :: n,m
            character (len=*), intent(inout), target :: a(n,m)
            integer (kind=IntKind), intent(in) :: msgtyp
            integer, intent(in) :: source
            integer (kind=IntKind) :: id
            integer (kind=IntKind) :: la
            if (source == MyPE_comm) then
!              l=min(slen,len(a(n,m)))
!              do k=1,m
!                 do j=1,n
!                    do i=1,l
!                       a(j,k)(i:i)=as2(i,j,k)
!                    enddo
!                 enddo
!              enddo
!              deallocate( as2 )
               id=1
               return
            else if (source > NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            else if (.not.isParentAlive() .and. source == NumPEs_comm) then
               call ErrorHandler('nbrecvMessage','invalid sopurce',   &
     &                           source)
            endif
            if (NumWaits==MaxWaits) then
               call printWaitInfo()
               call ErrorHandler('nbrecv_str2',                      &
     &                           'Too many messages are waiting!')
            endif
            NumWaits=NumWaits+1
            NumRecvs=NumRecvs+1
            do i=1,MaxWaits
               if (WaitType(i)=='NONE') then
                  id=i
                  exit
               endif
            enddo
            WaitType(id)='RECV'
            la=len(a(n,m))
#ifdef MPI
            if (source.ge.0) then
!              -------------------------------------------------------
               call MPI_irecv(a,len(a(1,1))*n*m,MPI_CHARACTER,source, &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=source
            else
!              -------------------------------------------------------
               call MPI_irecv(a,len(a(1,1))*n*m,MPI_CHARACTER, &
     &                        MPI_ANY_SOURCE, &
     &                        msgtyp,communicator,MsgID(id),info)
!              -------------------------------------------------------
               sourceNode=-1
            endif
#endif
         end function nbrecv_str2
!=====================================================================
!=====================================================================
!=====================================================================
         subroutine waitMessage(id)
            implicit none
            integer (kind=IntKind), intent(in) :: id
!
            if (.not.isParentAlive() .and. NumPEs_comm==1) then
               return
            else if (id<1 .or. id>MaxWaits) then
               call ErrorHandler('waitMessage', &
     &                           'invalid message id = ',id)
            else if (WaitType(id)=='RECV') then
               if (NumRecvs==0) then
                  call WarningHandler('waitMessage',                  &
     &                                'No message in waiting')
                  WaitType(id) = 'NONE'
                  return
               endif
               NumRecvs=NumRecvs-1
            else if (WaitType(id)=='SEND' .or. WaitType(id)=='SDPG') then
               if (NumSends==0) then
                  call WarningHandler('waitMessage',                  &
     &                                'No message in waiting')
                  WaitType(id) = 'NONE'
                  return
               endif
               NumSends=NumSends-1
            else
               call ErrorHandler('waitMessage', &
     &                           'invalid message id = ',id)
            endif
#ifdef MPI
            call MPI_wait(MsgID(id),status,info)
!
            if (WaitType(id)=='RECV') then
               sourceNode=status(MPI_SOURCE)
            endif
#endif
            MsgID(id)=0
            MsgSRC(id)=-1
            MsgBuffer(id)=0
            WaitType(id)='NONE'
            NumWaits=NumWaits-1
         end subroutine waitMessage
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
         subroutine printWaitInfo()
            implicit none
            integer (kind=IntKind) :: i
            write(6,'(3(a,i5,4x))')'MyPE = ',MyPE_comm,                   &
     &                             'MaxWaits = ',MaxWaits,           &
     &                             'NumWaits = ',NumWaits
            do i=1,NumWaits
               write(6,'(2(a,i5,4x),a,i10)')'MyPE = ',MyPE_comm,          &
     &                           'Index = ',i,'MsgID = ',MsgID(i)
               write(6,'(3(a,i5,4x))')'MyPE = ',MyPE_comm,                &
     &                           'Index = ',i,'MsgBuf = ',MsgBuffer(i)
               write(6,'(3(a,i5,4x))')'MyPE = ',MyPE_comm,                &
     &                           'Index = ',i,'MsgSRC = ',MsgSRC(i)
               write(6,'(2(a,i5,4x),a,a)')'MyPE = ',MyPE_comm,            &
     &                         'Index = ',i,'WaitType = ',WaitType(i)
            enddo
         end subroutine printWaitInfo
!=====================================================================
!
         function getNumWaits() result(n)
            implicit none
            integer (kind=IntKind) :: n
            n = NumWaits
         end function getNumWaits
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      subroutine openLocalMemory_s(a,n,accessID)
!=====================================================================
      implicit none
      integer (kind=IntKind), intent(in) :: n
      integer (kind=IntKind), intent(out) :: accessID
      integer (kind=IntKind) :: ierr
      integer (kind=MPI_ADDRESS_KIND) :: m
!
      character (len=*), intent(in) :: a(n)
!
#ifdef MPI
      m = len(a(1))*n*CHARACTER_BYTES
!     ----------------------------------------------------------------
      call MPI_win_create(a, m, CHARACTER_BYTES, MPI_INFO_NULL,       &
     &                    communicator, accessID, ierr)
!     ----------------------------------------------------------------
#else
      accessID = -1
#endif
!
      end subroutine openLocalMemory_s
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      subroutine openLocalMemory_s1(a,n1,n2,accessID)
!=====================================================================
      implicit none
      integer (kind=IntKind), intent(in) :: n1, n2
      integer (kind=IntKind), intent(out) :: accessID
      integer (kind=IntKind) :: ierr
      integer (kind=MPI_ADDRESS_KIND) :: m
!
      character (len=*), intent(in) :: a(n1,n2)
!
#ifdef MPI
      m = len(a(1,1))*n1*n2*CHARACTER_BYTES
!     ----------------------------------------------------------------
      call MPI_win_create(a, m, CHARACTER_BYTES, MPI_INFO_NULL,       &
     &                    communicator, accessID, ierr)
!     ----------------------------------------------------------------
#else
      accessID = -1
#endif
!
      end subroutine openLocalMemory_s1
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      subroutine openLocalMemory_i(a,n,accessID)
!=====================================================================
      implicit none
      integer (kind=IntKind), intent(in) :: n
      integer (kind=IntKind), intent(out) :: accessID
      integer (kind=IntKind) :: ierr
      integer (kind=MPI_ADDRESS_KIND) :: m
!
      integer (kind=IntKind), intent(in) :: a(n)
!
#ifdef MPI
      m = n*INTEGER_BYTES
!     ----------------------------------------------------------------
      call MPI_win_create(a, m, INTEGER_BYTES, MPI_INFO_NULL,         &
     &                    communicator, accessID, ierr)
!     ----------------------------------------------------------------
#else
      accessID = -1
#endif
!
      end subroutine openLocalMemory_i
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      subroutine openLocalMemory_i1(a,n1,n2,accessID)
!=====================================================================
      implicit none
      integer (kind=IntKind), intent(in) :: n1, n2
      integer (kind=IntKind), intent(out) :: accessID
      integer (kind=IntKind) :: ierr
      integer (kind=MPI_ADDRESS_KIND) :: m
!
      integer (kind=IntKind), intent(in) :: a(n1,n2)
!
#ifdef MPI
      m = n1*n2*INTEGER_BYTES
!     ----------------------------------------------------------------
      call MPI_win_create(a, m, INTEGER_BYTES, MPI_INFO_NULL,         &
     &                    communicator, accessID, ierr)
!     ----------------------------------------------------------------
#else
      accessID = -1
#endif
!
      end subroutine openLocalMemory_i1
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      subroutine openLocalMemory_r(a,n,accessID)
!=====================================================================
      implicit none
      integer (kind=IntKind), intent(in) :: n
      integer (kind=IntKind), intent(out) :: accessID
      integer (kind=IntKind) :: ierr
      integer (kind=MPI_ADDRESS_KIND) :: m
!
      real (kind=RealKind), intent(in) :: a(n)
!
#ifdef MPI
      m = n*REAL_BYTES
!     ----------------------------------------------------------------
      call MPI_win_create(a, m, REAL_BYTES, MPI_INFO_NULL,            &
     &                    communicator, accessID, ierr)
!     ----------------------------------------------------------------
#else
      accessID = -1
#endif
!
      end subroutine openLocalMemory_r
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      subroutine openLocalMemory_r1(a,n1,n2,accessID)
!=====================================================================
      implicit none
      integer (kind=IntKind), intent(in) :: n1, n2
      integer (kind=IntKind), intent(out) :: accessID
      integer (kind=IntKind) :: ierr
      integer (kind=MPI_ADDRESS_KIND) :: m
!
      real (kind=RealKind), intent(in) :: a(n1,n2)
!
#ifdef MPI
      m = n1*n2*REAL_BYTES
!     ----------------------------------------------------------------
      call MPI_win_create(a, m, REAL_BYTES, MPI_INFO_NULL,            &
     &                    communicator, accessID, ierr)
!     ----------------------------------------------------------------
#else
      accessID = -1
#endif
!
      end subroutine openLocalMemory_r1
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      subroutine openLocalMemory_c(a,n,accessID)
!=====================================================================
      implicit none
      integer (kind=IntKind), intent(in) :: n
      integer (kind=IntKind), intent(out) :: accessID
      integer (kind=IntKind) :: ierr
      integer (kind=MPI_ADDRESS_KIND) :: m
!
      complex (kind=CmplxKind), intent(in) :: a(n)
!
#ifdef MPI
      m = n*COMPLEX_BYTES
!     ----------------------------------------------------------------
      call MPI_win_create(a, m, COMPLEX_BYTES, MPI_INFO_NULL,         &
     &                    communicator, accessID, ierr)
!     ----------------------------------------------------------------
#else
      accessID = -1
#endif
!
      end subroutine openLocalMemory_c
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      subroutine openLocalMemory_c1(a,n1,n2,accessID)
!=====================================================================
      implicit none
      integer (kind=IntKind), intent(in) :: n1, n2
      integer (kind=IntKind), intent(out) :: accessID
      integer (kind=IntKind) :: ierr
      integer (kind=MPI_ADDRESS_KIND) :: m
!
      complex (kind=CmplxKind), intent(in) :: a(n1,n2)
!
#ifdef MPI
      m = n1*n2*COMPLEX_BYTES
!     ----------------------------------------------------------------
      call MPI_win_create(a, m, COMPLEX_BYTES, MPI_INFO_NULL,         &
     &                    communicator, accessID, ierr)
!     ----------------------------------------------------------------
#else
      accessID = -1
#endif
!
      end subroutine openLocalMemory_c1
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      subroutine readRemoteMemory_s(a,n,proc,accessID,index)
!=====================================================================
      implicit none
      integer (kind=IntKind), intent(in) :: n
      integer (kind=IntKind), intent(in) :: accessID
      integer (kind=IntKind), intent(in) :: proc
      integer (kind=IntKind), intent(in) :: index
      integer (kind=IntKind) :: ierr, m
      integer (kind=MPI_ADDRESS_KIND) :: offset
!
      character (len=*), intent(out) :: a(n)
!
#ifdef MPI
!     ----------------------------------------------------------------
#ifndef No_MPI_LOCK
      call MPI_win_lock(MPI_LOCK_SHARED, proc, 0, accessID, ierr)
#endif
      offset = index-1
      m = n*len(a(1))
      call MPI_get(a, m, MPI_CHARACTER, proc, offset,                &
     &             m, MPI_CHARACTER, accessID, ierr)
#ifndef No_MPI_LOCK
      call MPI_win_unlock(proc, accessID, ierr)
#endif
!     ----------------------------------------------------------------
#else
      a(1:n) = ' '
#endif
!
      end subroutine readRemoteMemory_s
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      subroutine readRemoteMemory_i(a,n,proc,accessID,index)
!=====================================================================
      implicit none
      integer (kind=IntKind), intent(in) :: n
      integer (kind=IntKind), intent(in) :: accessID
      integer (kind=IntKind), intent(in) :: proc
      integer (kind=IntKind), intent(in) :: index
      integer (kind=IntKind) :: ierr
      integer (kind=MPI_ADDRESS_KIND) :: offset
!
      integer (kind=IntKind), intent(out) :: a(n)
!
#ifdef MPI
!     ----------------------------------------------------------------
#ifndef No_MPI_LOCK
      call MPI_win_lock(MPI_LOCK_SHARED, proc, 0, accessID, ierr)
#endif
      offset = index-1
      call MPI_get(a, n, MPI_INTEGER, proc, offset,                &
     &             n, MPI_INTEGER, accessID, ierr)
#ifndef No_MPI_LOCK
      call MPI_win_unlock(proc, accessID, ierr)
#endif
!     ----------------------------------------------------------------
#else
      a(1:n) = 0
#endif
!
      end subroutine readRemoteMemory_i
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      subroutine readRemoteMemory_r(a,n,proc,accessID,index)
!=====================================================================
      implicit none
      integer (kind=IntKind), intent(in) :: n
      integer (kind=IntKind), intent(in) :: accessID
      integer (kind=IntKind), intent(in) :: proc
      integer (kind=IntKind), intent(in) :: index
      integer (kind=IntKind) :: ierr
      integer (kind=MPI_ADDRESS_KIND) :: offset
!
      real (kind=RealKind), intent(out) :: a(n)
!
#ifdef MPI
!     ----------------------------------------------------------------
#ifndef No_MPI_LOCK
      call MPI_win_lock(MPI_LOCK_SHARED, proc, 0, accessID, ierr)
#endif
      offset = index-1
      call MPI_get(a, n, MPI_DOUBLE_PRECISION, proc, offset,         &
     &             n, MPI_DOUBLE_PRECISION, accessID, ierr)
#ifndef No_MPI_LOCK
      call MPI_win_unlock(proc, accessID, ierr)
#endif
!     ----------------------------------------------------------------
#else
      a(1:n) = 0.0d0
#endif
!
      end subroutine readRemoteMemory_r
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      subroutine readRemoteMemory_c(a,n,proc,accessID,index)
!=====================================================================
      implicit none
      integer (kind=IntKind), intent(in) :: n
      integer (kind=IntKind), intent(in) :: accessID
      integer (kind=IntKind), intent(in) :: proc
      integer (kind=IntKind), intent(in) :: index
      integer (kind=IntKind) :: ierr
      integer (kind=MPI_ADDRESS_KIND) :: offset
!
      complex (kind=CmplxKind), intent(out) :: a(n)
!
#ifdef MPI
!     ----------------------------------------------------------------
#ifndef No_MPI_LOCK
      call MPI_win_lock(MPI_LOCK_SHARED, proc, 0, accessID, ierr)
#endif
      offset = index-1
      call MPI_get(a, n, MPI_DOUBLE_COMPLEX, proc, offset,           &
     &             n, MPI_DOUBLE_COMPLEX, accessID, ierr)
#ifndef No_MPI_LOCK
      call MPI_win_unlock(proc, accessID, ierr)
#endif
!     ----------------------------------------------------------------
#else
      a(1:n) = cmplx(0.0d0,0.0d0,kind=CmplxKind)
#endif
!
      end subroutine readRemoteMemory_c
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      subroutine writeRemoteMemory_s(a,n,proc,accessID,index)
!=====================================================================
      implicit none
      integer (kind=IntKind), intent(in) :: n
      integer (kind=IntKind), intent(in) :: accessID
      integer (kind=IntKind), intent(in) :: proc
      integer (kind=IntKind), intent(in) :: index
      integer (kind=IntKind) :: ierr, m
      integer (kind=MPI_ADDRESS_KIND) :: offset
!
      character (len=*), intent(in) :: a(n)
!
#ifdef MPI
!     ----------------------------------------------------------------
#ifndef No_MPI_LOCK
      call MPI_win_lock(MPI_LOCK_SHARED, proc, 0, accessID, ierr)
#endif
      offset = index-1
      m = n*len(a(1))
      call MPI_put(a, m, MPI_CHARACTER, proc, offset,                &
     &             m, MPI_CHARACTER, accessID, ierr)
#ifndef No_MPI_LOCK
      call MPI_win_unlock(proc, accessID, ierr)
#endif
!     ----------------------------------------------------------------
#endif
!
      end subroutine writeRemoteMemory_s
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      subroutine writeRemoteMemory_i(a,n,proc,accessID,index)
!=====================================================================
      implicit none
      integer (kind=IntKind), intent(in) :: n
      integer (kind=IntKind), intent(in) :: accessID
      integer (kind=IntKind), intent(in) :: proc
      integer (kind=IntKind), intent(in) :: index
      integer (kind=IntKind) :: ierr
      integer (kind=MPI_ADDRESS_KIND) :: offset
!
      integer (kind=IntKind), intent(in) :: a(n)
!
#ifdef MPI
!     ----------------------------------------------------------------
#ifndef No_MPI_LOCK
      call MPI_win_lock(MPI_LOCK_SHARED, proc, 0, accessID, ierr)
#endif
      offset = index-1
      call MPI_put(a, n, MPI_INTEGER, proc, offset,                  &
     &             n, MPI_INTEGER, accessID, ierr)
#ifndef No_MPI_LOCK
      call MPI_win_unlock(proc, accessID, ierr)
#endif
!     ----------------------------------------------------------------
#endif
!
      end subroutine writeRemoteMemory_i
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      subroutine writeRemoteMemory_r(a,n,proc,accessID,index)
!=====================================================================
      implicit none
      integer (kind=IntKind), intent(in) :: n
      integer (kind=IntKind), intent(in) :: accessID
      integer (kind=IntKind), intent(in) :: proc
      integer (kind=IntKind), intent(in) :: index
      integer (kind=IntKind) :: ierr
      integer (kind=MPI_ADDRESS_KIND) :: offset
!
      real (kind=RealKind), intent(in) :: a(n)
!
#ifdef MPI
!     ----------------------------------------------------------------
#ifndef No_MPI_LOCK
      call MPI_win_lock(MPI_LOCK_SHARED, proc, 0, accessID, ierr)
#endif
      offset = index-1
      call MPI_put(a, n, MPI_DOUBLE_PRECISION, proc, offset,         &
     &             n, MPI_DOUBLE_PRECISION, accessID, ierr)
#ifndef No_MPI_LOCK
      call MPI_win_unlock(proc, accessID, ierr)
#endif
!     ----------------------------------------------------------------
#endif
!
      end subroutine writeRemoteMemory_r
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      subroutine writeRemoteMemory_c(a,n,proc,accessID,index)
!=====================================================================
      implicit none
      integer (kind=IntKind), intent(in) :: n
      integer (kind=IntKind), intent(in) :: accessID
      integer (kind=IntKind), intent(in) :: proc
      integer (kind=IntKind), intent(in) :: index
      integer (kind=IntKind) :: ierr
      integer (kind=MPI_ADDRESS_KIND) :: offset
!
      complex (kind=CmplxKind), intent(in) :: a(n)
!
#ifdef MPI
!     ----------------------------------------------------------------
#ifndef No_MPI_LOCK
      call MPI_win_lock(MPI_LOCK_SHARED, proc, 0, accessID, ierr)
#endif
      offset = index-1
      call MPI_put(a, n, MPI_DOUBLE_COMPLEX, proc, offset,           &
     &             n, MPI_DOUBLE_COMPLEX, accessID, ierr)
#ifndef No_MPI_LOCK
      call MPI_win_unlock(proc, accessID, ierr)
#endif
!     ----------------------------------------------------------------
#endif
!
      end subroutine writeRemoteMemory_c
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      subroutine closeLocalMemory(AccessID)
!=====================================================================
      implicit none
      integer (kind=IntKind), intent(in) :: AccessID
      integer (kind=IntKind) :: ierr
!
#ifdef MPI
!     ----------------------------------------------------------------
      call MPI_win_free(AccessID, ierr)
!     ----------------------------------------------------------------
#endif
      end subroutine closeLocalMemory
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      subroutine syncLocalMemory(AccessID)
!=====================================================================
      implicit none
      integer (kind=IntKind), intent(in) :: AccessID
      integer (kind=IntKind) :: ierr
!
#ifdef MPI
#ifdef No_MPI_LOCK
!     ----------------------------------------------------------------
      call MPI_WIN_FENCE(0, accessID, ierr)
!     ----------------------------------------------------------------
#endif
#endif
      end subroutine syncLocalMemory
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      function getCommunicator() result(comm)
      implicit none
!
      integer (kind=IntKind) :: comm
!
      comm = communicator
!
      end function getCommunicator
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      subroutine setCommunicator(comm,p,n,sync)
      implicit none
!
      integer (kind=IntKind), intent(in) :: comm, p, n
      logical, intent(in), optional :: sync
!
      if (p >= n .or. p < 0 .or. n < 0) then
         call ErrorHandler('setCommunicator',                         &
              'Invalid MyPE ond/or NumPEs input for the new comm',    &
                           p, n)
      endif
!
      communicator = comm
      MyPE_comm = p
      NumPEs_comm = n
!
#ifdef MPI
      if (present(sync)) then
         if (sync) then
            call MPI_barrier(communicator,info)
         endif
      endif
#endif
!
      end subroutine setCommunicator
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      subroutine resetCommunicator(sync)
      implicit none
!
      logical, intent(in), optional :: sync
!
#ifdef MPI
      if (present(sync)) then
         if (sync) then
            call MPI_barrier(communicator,info)
         endif
      endif
#endif
!
      communicator = MPI_COMM_WORLD
      MyPE_comm = MyPE
      NumPEs_comm = NumPEs
!
      end subroutine resetCommunicator
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      function isThereMessage(msgtyp,source,blocking) result(y)
      implicit none
      integer (kind=IntKind), intent(in) :: msgtyp
      integer (kind=IntKind), intent(in) :: source
      logical, optional, intent(in) :: blocking
      logical :: y
!
#ifdef MPI
      if (present(blocking)) then
         if (blocking) then
!           ----------------------------------------------------------
            call MPI_probe(source,msgtyp,communicator,status,info)
!           ----------------------------------------------------------
            y = .true.
            sourceNode=status(MPI_SOURCE)
            return
         endif
      endif
!
      y = .false.
!     ----------------------------------------------------------------
      call MPI_iprobe(source,msgtyp,communicator,y,status,info)
!     ----------------------------------------------------------------
      if (y) then
         sourceNode=status(MPI_SOURCE)
      endif
#else
      y = .false.
#endif
!
      end function isThereMessage
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      function getMaxWaits() result(n)
      implicit none
!
      integer (kind=IntKind) :: n
!
      n = MaxWaits
!
      end function getMaxWaits
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      subroutine setMaxWaits(n)
      implicit none
!
      integer (kind=IntKind), intent(in) :: n
      integer (kind=IntKind), allocatable :: MsgIDTmp(:)
      integer (kind=IntKind), allocatable :: MsgBufferTmp(:)
      integer (kind=IntKind), allocatable :: MsgSRCTmp(:)
      character (len=4), allocatable :: WaitTypeTmp(:)
!
      if (n <= MaxWaits) then
!        call WarningHandler('setMaxWaits','MaxWaits is unchanged for n < MaxWaits', &
!                            n, MaxWaits)
         return
      endif
!
      allocate(MsgIDTmp(MaxWaits), MsgBufferTmp(MaxWaits),            &
               MsgSRCTmp(MaxWaits), WaitTypeTmp(MaxWaits))
      MsgIDTmp(1:MaxWaits) = MsgID(1:MaxWaits)
      MsgBufferTmp(1:MaxWaits) = MsgBuffer(1:MaxWaits)
      MsgSRCTmp(1:MaxWaits) = MsgSRC(1:MaxWaits)
      WaitTypeTmp(1:MaxWaits) = WaitType(1:MaxWaits)
!
      deallocate(MsgID, MsgBuffer, MsgSRC, WaitType)
      allocate(MsgID(n), MsgBuffer(n), MsgSRC(n), WaitType(n))
!
      MsgID(MaxWaits+1:n)=0
      MsgBuffer(MaxWaits+1:n)=0
      MsgSRC(MaxWaits+1:n)=-1
      WaitType(MaxWaits+1:n)='NONE'
      MsgID(1:MaxWaits) = MsgIDTmp(1:MaxWaits)
      MsgBuffer(1:MaxWaits) = MsgBufferTmp(1:MaxWaits)
      MsgSRC(1:MaxWaits) = MsgSRCTmp(1:MaxWaits)
      WaitType(1:MaxWaits) = WaitTypeTmp(1:MaxWaits)
!
      MaxWaits = n
!
      deallocate(MsgIDTmp, MsgBufferTmp, MsgSRCTmp, WaitTypeTmp)
!
      end subroutine setMaxWaits
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      subroutine startRoundTurn(pe)
      implicit none
      integer (kind=IntKind), intent(in), optional :: pe
      integer (kind=IntKind) :: msg, source, wait_id
!
      if (present(pe)) then
         if (pe < 0 .or. pe > NumPEs-1) then
            call ErrorHandler('startRoundTurn','starting PE is out of range',pe)
         endif
         first_robin = pe
      else
         first_robin = 0
      endif
!
      if (MyPE /= first_robin) then
         source = MyPE - 1
         if (source < 0) then
            source = NumPEs - 1
         endif
#ifdef MPI
!        -------------------------------------------------------------
         call MPI_irecv(msg,1,MPI_INTEGER,source,                     &
     &                  robin_message_type,communicator,wait_id,info)
!        -------------------------------------------------------------
         call MPI_wait(wait_id,status,info)
!        -------------------------------------------------------------
#endif
      endif
!
      end subroutine startRoundTurn
!=====================================================================
!
!*********************************************************************
!
!=====================================================================
      subroutine finishMyTurn()
      implicit none
      integer (kind=IntKind) :: msg, next, wait_id
!
      next = MyPE + 1
      if (next >= NumPEs) then
         next = 0
      endif
#ifdef MPI
      if (next /= first_robin) then
         msg = MyPE
!        -------------------------------------------------------------
         call MPI_isend(msg,1,MPI_INTEGER,next,                       &
     &                  robin_message_type,communicator,wait_id,info)
         call MPI_wait(wait_id,status,info)
!        -------------------------------------------------------------
      endif
#endif
!
      end subroutine finishMyTurn
!=====================================================================
      end module MPPModule
