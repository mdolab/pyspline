/*
 * NOTICE and LICENSE for Tecplot Input/Output Library (TecIO) - OpenFOAM
 *
 * Copyright (C) 1988-2009 Tecplot, Inc.  All rights reserved worldwide.
 *
 * Tecplot hereby grants OpenCFD limited authority to distribute without
 * alteration the source code to the Tecplot Input/Output library, known 
 * as TecIO, as part of its distribution of OpenFOAM and the 
 * OpenFOAM_to_Tecplot converter.  Users of this converter are also hereby
 * granted access to the TecIO source code, and may redistribute it for the
 * purpose of maintaining the converter.  However, no authority is granted
 * to alter the TecIO source code in any form or manner.
 *
 * This limited grant of distribution does not supersede Tecplot, Inc.'s 
 * copyright in TecIO.  Contact Tecplot, Inc. for further information.
 * 
 * Tecplot, Inc.
 * 3535 Factoria Blvd, Ste. 550
 * Bellevue, WA 98006, USA
 * Phone: +1 425 653 1200
 * http://www.tecplot.com/
 *
 */
#include "stdafx.h"
#include "MASTER.h"
#define TECPLOTENGINEMODULE

/*
 *****************************************************************
 *****************************************************************
 *******                                                  ********
 ****** Copyright (C) 1988-2008 Tecplot, Inc.               ********
 *******       All Rights Reserved.                       ********
 *******                                                  ********
 *****************************************************************
 *****************************************************************
 */

#define ARRLISTMODULE
#include "GLOBAL.h"
#include "TASSERT.h"
#include "ALLOC.h"
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
#include "ARRLIST.h"

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

/*
 * ABSTRACT:
 *
 * This general purpose list uses an array implementation. The high use member
 * functions have macro covers to make the implementation both efficient with
 * respect to speed without compromising the internal representation. The
 * internal array is allocated to fit the requested type's item size. Most
 * intrinsic 'C' and Tecplot types are available.
 */


/**
 * Copies the private array items from the specified source to the target
 * location. The buffers may overlap.
 *
 * note
 *     Originally this function was a macro that called memmove
 *     directly:
 *
 *         #define CopyArrayItems(TargetArray, TargetOffset, \
 *                                SourceArray, SourceOffset, \
 *                                Count, ItemSize) \
 *                     (memmove(&((TargetArray)[(TargetOffset)*ItemSize]), \
 *                              &((SourceArray)[(SourceOffset)*ItemSize]), \
 *                              Count*ItemSize))
 *
 * This however proved troublesome as some machines replaced the memmove
 * with a call to memcpy in the linker. The memcpy function does not support
 * overlapping moves so I could not use it. This function should be just
 * about as fast however so it is no big deal.
 *
 * param TargetArray
 *     Base address of the target array to receive the items.
 * param TargetOffset
 *     Target offset of the first item.
 * param SourceArray
 *     Base address of the source array supplying the items.
 * param SourceOffset
 *     Source offset of the first item.
 * param Count
 *     Number of items to copy.
 * param ItemSize
 *     Item size in bytes.
 */
static void CopyArrayItems(char        *TargetArray,
                           LgIndex_t    TargetOffset,
                           char        *SourceArray,
                           LgIndex_t    SourceOffset,
                           LgIndex_t    Count,
                           SmInteger_t  ItemSize)
{
    REQUIRE(VALID_REF(TargetArray));
    REQUIRE(TargetOffset >= 0);
    REQUIRE(VALID_REF(SourceArray));
    REQUIRE(SourceOffset >= 0);
    REQUIRE(&TargetArray[TargetOffset] != &SourceArray[SourceOffset]);
    REQUIRE(Count >= 1);
    REQUIRE(1 <= ItemSize && ItemSize <= (SmInteger_t)sizeof(ArrayListItem_u));

    void* TargetPtr = &TargetArray[TargetOffset * ItemSize];
    void* SourcePtr = &SourceArray[SourceOffset * ItemSize];
    memmove(TargetPtr, SourcePtr, ((size_t)Count) * ItemSize);
}


/**
 * Adjusts the capacity request as necessary to minimize memory reallocations
 * for large lists. Unless the request exceeds the maximum the adjusted
 * capacity will be at least as big as requested however it may be larger if it
 * is determined that the space requirement is growing faster. If the maximum
 * is exceeded zero should be returned.
 *
 * param ArrayList
 *     Array list requesting the change in capacity.
 * param CurrentCapacity
 *     Current capacity of the array list.
 * param RequestedCapacity
 *     Capacity request or zero for default size.
 * param ClientData
 *     Any client data needed for the adjustment.
 *
 * return
 *     Adjusted capacity that is at least as large as the request or zero if
 *     unable to satisfy the requested capacity.
 */
static LgIndex_t AdjustCapacityRequest(ArrayList_pa ArrayList,
                                       LgIndex_t    CurrentCapacity,
                                       LgIndex_t    RequestedCapacity,
                                       ArbParam_t   ClientData)
{
    LgIndex_t Result;

    REQUIRE(ArrayListIsValid(ArrayList));
    REQUIRE((RequestedCapacity == 0 && CurrentCapacity == 0) ||
            RequestedCapacity > CurrentCapacity);

    if (RequestedCapacity != 0 && CurrentCapacity == 0)
    {
        /* first allocation; assume the request is the desired capacityy */
        Result = RequestedCapacity;
    }
    else
    {
        const LgIndex_t DEFAULT_CAPACITY = 32;
        LgIndex_t       BlockSize = MAX(DEFAULT_CAPACITY, CurrentCapacity / 2);
        if (RequestedCapacity == 0)
            Result = DEFAULT_CAPACITY;
        else
            Result = ((RequestedCapacity - 1) / BlockSize + 1) * BlockSize;
    }

    ENSURE(Result == 0 || Result >= RequestedCapacity);
    return Result;
}


/**
 * Gets the size of an individual element.
 *
 * param Type
 *     Array list element type.
 *
 * return
 *     Element size corresponding to the type.
 */
static SmInteger_t GetElementSize(ArrayListType_e Type)
{
    SmInteger_t Result;

    REQUIRE(VALID_ENUM(Type, ArrayListType_e));

    switch (Type)
    {
        case ArrayListType_UnsignedChar:
            Result = (SmInteger_t)sizeof(unsigned char);
            break;
        case ArrayListType_UnsignedShort:
            Result = (SmInteger_t)sizeof(unsigned short);
            break;
        case ArrayListType_UnsignedInt:
            Result = (SmInteger_t)sizeof(unsigned int);
            break;
        case ArrayListType_UnsignedLong:
            Result = (SmInteger_t)sizeof(unsigned long);
            break;
        case ArrayListType_Int64:
            Result = (SmInteger_t)sizeof(Int64_t);
            break;
        case ArrayListType_Char:
            Result = (SmInteger_t)sizeof(char);
            break;
        case ArrayListType_Short:
            Result = (SmInteger_t)sizeof(short);
            break;
        case ArrayListType_Int:
            Result = (SmInteger_t)sizeof(int);
            break;
        case ArrayListType_Long:
            Result = (SmInteger_t)sizeof(long);
            break;
        case ArrayListType_Float:
            Result = (SmInteger_t)sizeof(float);
            break;
        case ArrayListType_Double:
            Result = (SmInteger_t)sizeof(double);
            break;
        case ArrayListType_LgIndex:
            Result = (SmInteger_t)sizeof(LgIndex_t);
            break;
        case ArrayListType_EntIndex:
            Result = (SmInteger_t)sizeof(EntIndex_t);
            break;
        case ArrayListType_SmInteger:
            Result = (SmInteger_t)sizeof(SmInteger_t);
            break;
        case ArrayListType_Boolean:
            Result = (SmInteger_t)sizeof(Boolean_t);
            break;
        case ArrayListType_ArbParam:
            Result = (SmInteger_t)sizeof(ArbParam_t);
            break;

        case ArrayListType_UnsignedCharPtr:
            Result = (SmInteger_t)sizeof(unsigned char *);
            break;
        case ArrayListType_UnsignedShortPtr:
            Result = (SmInteger_t)sizeof(unsigned short *);
            break;
        case ArrayListType_UnsignedIntPtr:
            Result = (SmInteger_t)sizeof(unsigned int *);
            break;
        case ArrayListType_UnsignedLongPtr:
            Result = (SmInteger_t)sizeof(unsigned long *);
            break;
        case ArrayListType_Int64Ptr:
            Result = (SmInteger_t)sizeof(Int64_t *);
            break;
        case ArrayListType_CharPtr:
            Result = (SmInteger_t)sizeof(char *);
            break;
        case ArrayListType_ShortPtr:
            Result = (SmInteger_t)sizeof(short *);
            break;
        case ArrayListType_IntPtr:
            Result = (SmInteger_t)sizeof(int *);
            break;
        case ArrayListType_LongPtr:
            Result = (SmInteger_t)sizeof(long *);
            break;
        case ArrayListType_FloatPtr:
            Result = (SmInteger_t)sizeof(float *);
            break;
        case ArrayListType_DoublePtr:
            Result = (SmInteger_t)sizeof(double *);
            break;
        case ArrayListType_LgIndexPtr:
            Result = (SmInteger_t)sizeof(LgIndex_t *);
            break;
        case ArrayListType_EntIndexPtr:
            Result = (SmInteger_t)sizeof(EntIndex_t *);
            break;
        case ArrayListType_SmIntegerPtr:
            Result = (SmInteger_t)sizeof(SmInteger_t *);
            break;
        case ArrayListType_BooleanPtr:
            Result = (SmInteger_t)sizeof(Boolean_t *);
            break;
        case ArrayListType_ArbParamPtr:
            Result = (SmInteger_t)sizeof(ArbParam_t *);
            break;

        case ArrayListType_VoidPtr:
            Result = (SmInteger_t)sizeof(void *);
            break;
        case ArrayListType_FunctionPtr:
            Result = (SmInteger_t)sizeof(void (*)());
            break;
        case ArrayListType_Any: /* allows a mixed bag of items */
            Result = (SmInteger_t)sizeof(ArrayListItem_u);
            break;

        default:
            Result = 0; /* make some compilers happy. */
            CHECK(FALSE);
            break;
    }

    ENSURE(1 <= Result && Result <= (SmInteger_t)sizeof(ArrayListItem_u));
    return Result;
}


/**
 * Calls the item destructor for each item specified.
 *
 * param ArrayList
 *     Array list needing its items destroyed.
 * param ItemOffset
 *     Offset to the first item to destroy in the list.
 * param ItemSize
 *     Size of each array list item.
 * param Count
 *     Number of items to destroy.
 * param ItemDestructor
 *     Function called for each array list item.
 * param CientData
 *     Any client data needed for the destructor.
 */
static void DestroyItems(ArrayList_pa               ArrayList,
                         LgIndex_t                  ItemOffset,
                         SmInteger_t                ItemSize,
                         LgIndex_t                  Count,
                         ArrayListItemDestructor_pf ItemDestructor,
                         ArbParam_t                 ClientData)
{
    LgIndex_t Index;

    REQUIRE(ArrayListIsValid(ArrayList));
    REQUIRE(0 <= ItemOffset && ItemOffset <= ArrayList->Count - 1);
    REQUIRE(1 <= Count && ItemOffset + Count <= ArrayList->Count);
    REQUIRE(VALID_FN_REF(ItemDestructor));

    for (Index = 0;
         Index < Count;
         Index++)
    {
        LgIndex_t ItemIndex = (Index + ItemOffset) * ItemSize;
#if !defined NO_ASSERTS
        Boolean_t DoContinue = ItemDestructor((void *) & ArrayList->Array[ItemIndex], ClientData);
        CHECK(DoContinue); /* this is a requirement of ArrayListItemDestructor_pf */
#else
        ItemDestructor((void *)&ArrayList->Array[ItemIndex], ClientData);
#endif
    }
}


/**
 * Calls the item duplicator for each item specified.
 *
 * param TargetArray
 *     Target array needing its items duplicated.
 * param TargetItemOffset
 *     Target offset to the first duplicated item.
 * param SourceArray
 *     Source array needing its items duplicated.
 * param SourceItemOffset
 *     Source offset to the first item to duplicate in the list.
 * param ItemSize
 *     Size of each array list item.
 * param Count
 *     Number of items to duplicate.
 * param ItemDuplicator
 *     Function called for each array list item.
 * param CientData
 *     Any client data needed for the destructor.
 *
 * return
 *     TRUE if the duplication was a success
 *     FALSE otherwise
 */
static Boolean_t DuplicateItems(char                       *TargetArray,
                                LgIndex_t                   TargetItemOffset,
                                char                       *SourceArray,
                                LgIndex_t                   SourceItemOffset,
                                SmInteger_t                 ItemSize,
                                LgIndex_t                   Count,
                                ArrayListItemDuplicator_pf  ItemDuplicator,
                                ArbParam_t                  ClientData)
{
    Boolean_t IsOk = TRUE;
    LgIndex_t Index;

    REQUIRE(VALID_REF(TargetArray));
    REQUIRE(TargetItemOffset >= 0);
    REQUIRE(VALID_REF(SourceArray));
    REQUIRE(SourceItemOffset >= 0);
    REQUIRE(1 <= ItemSize &&
            ItemSize <= (SmInteger_t)sizeof(ArrayListItem_u));
    REQUIRE(Count >= 1);
    REQUIRE(VALID_FN_REF(ItemDuplicator));

    for (Index = 0;
         Index < Count && IsOk;
         Index++)
    {
        LgIndex_t TargetItemIndex = (Index + TargetItemOffset) * ItemSize;
        LgIndex_t SourceItemIndex = (Index + SourceItemOffset) * ItemSize;
        IsOk = ItemDuplicator((void *) & TargetArray[TargetItemIndex],
                              (void *) & SourceArray[SourceItemIndex],
                              ClientData);
    }

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}


/**
 * Determine if the list handle is sane.
 *
 * param ArrayList
 *     Array list in question.
 *
 * return
 *     TRUE if the array list is valid, otherwise FALSE.
 */
Boolean_t ArrayListIsValid(ArrayList_pa ArrayList)
{
    Boolean_t IsValid;

    /* this just makes sure that the NULL item global was initialized */
    INVARIANT(ArrayListNullItem.Double == 0.0);

    IsValid = (VALID_REF(ArrayList) &&
               VALID_ENUM(ArrayList->Type, ArrayListType_e) &&
               (1 <= ArrayList->ItemSize &&
                ArrayList->ItemSize <= (SmInteger_t)sizeof(ArrayListItem_u)) &&
               (0 <= ArrayList->Count &&
                ArrayList->Count <= ArrayList->Capacity));

    ENSURE(VALID_BOOLEAN(IsValid));
    return IsValid;
}


/**
 * Gets the specified array list's type.
 *
 * param ArrayList
 *     Array list of which the type is desired.
 *
 * return
 *     Array list type.
 */
ArrayListType_e ArrayListGetType(ArrayList_pa ArrayList)
{
    ArrayListType_e Result;

    REQUIRE(ArrayListIsValid(ArrayList));

    Result = ArrayList->Type;

    ENSURE(VALID_ENUM(Result, ArrayListType_e));
    return Result;
}


/**
 * Enlarge the list capacity to accommodate, at a minimum, the requested
 * capacity (number of items). The enlarged section is initialized with zeros.
 * This function can be called by clients who want to ensure that calls to
 * the ArrayListSetXxx class of functions will never fail for offsets within
 * the RequestedCapacity.
 *
 * param ArrayList
 *     Current capacity used as a helpful hint for the adjustment algorythm.
 * param RequestedCapacity
 *     Capacity (number ot items) request or zero for default size.
 *
 * return
 *     TRUE if the list could be enlarged (or was large enough),
 *     otherwise FALSE.
 */
Boolean_t ArrayListEnlargeCapacity(ArrayList_pa ArrayList,
                                   LgIndex_t    RequestedCapacity)
{
    Boolean_t IsOk;
    LgIndex_t AdjustedCapacity;
    char      *EnlargedArray;

    REQUIRE(ArrayListIsValid(ArrayList));
    REQUIRE(IMPLICATION(RequestedCapacity == 0, ArrayList->Capacity == 0));

    if (RequestedCapacity == 0 ||
        RequestedCapacity > ArrayList->Capacity)
    {
        AdjustedCapacity =
            ArrayList->CapacityRequestAdjuster(ArrayList,
                                               ArrayList->Capacity,
                                               RequestedCapacity,
                                               ArrayList->CapacityRequestAdjusterClientData);
        CHECK(AdjustedCapacity == 0 ||
              AdjustedCapacity >= RequestedCapacity);

        IsOk = (AdjustedCapacity != 0); /* ...were we able to meet the request? */
        if (IsOk)
        {
            EnlargedArray = ALLOC_ARRAY(AdjustedCapacity * ArrayList->ItemSize,
                                        char, "array list");
            if (EnlargedArray == NULL)
            {
                /* try again with minimum capacity request */
                if (RequestedCapacity != 0)
                    AdjustedCapacity = RequestedCapacity;
                else
                    AdjustedCapacity = 1; /* can't get smaller than this */
                EnlargedArray = ALLOC_ARRAY(AdjustedCapacity * ArrayList->ItemSize,
                                            char, "array list");
            }
            IsOk = (EnlargedArray != NULL);
        }

        if (IsOk)
        {
            /*
             * Initialize the expanded section of the array with zeros. This default
             * value of zero is necessary for many other array list operations.
             */
            memset(&EnlargedArray[ArrayList->Count*ArrayList->ItemSize], 0,
                   (AdjustedCapacity - ArrayList->Count)*ArrayList->ItemSize);

            if (ArrayList->Array != NULL)
            {
                if (ArrayList->Count != 0)
                    CopyArrayItems(EnlargedArray, 0,
                                   ArrayList->Array, 0,
                                   ArrayList->Count,
                                   ArrayList->ItemSize);
                FREE_ARRAY(ArrayList->Array, "array list");
            }

            ArrayList->Array    = EnlargedArray;
            ArrayList->Capacity = AdjustedCapacity;
        }
    }
    else
    {
        IsOk = TRUE;
    }

    ENSURE(ArrayListIsValid(ArrayList));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}


/**
 * Allocates an array list handle with the estimated capacity
 * or a suitable default if an estimate is unavailable.
 *
 * param EstimatedCapacity
 *     Clients best guess at the estimated capacity (number of items) needed.
 *     If an estimate is not available zero the zero should be used to get the
 *     default capacity.
 * param Type
 *     Type of array list being allocated.
 * param CapacityRequestAdjuster
 *     Function to use to adjust any capacity change requests or
 *     NULL if the default adjuster is good enough.
 * param CapacityRequestAdjusterClientData
 *     Any client data needed for the capacity adjustment.
 *
 * return
 *     Array list handle if sufficient memory was available,
 *     otherwise a handle to NULL.
 */
ArrayList_pa ArrayListAlloc(LgIndex_t                           EstimatedCapacity,
                            ArrayListType_e                     Type,
                            ArrayListCapacityRequestAdjuster_pf CapacityRequestAdjuster,
                            ArbParam_t                          CapacityRequestAdjusterClientData)
{
    ArrayList_pa Result;

    REQUIRE(EstimatedCapacity >= 0);
    REQUIRE(VALID_ENUM(Type, ArrayListType_e));

    Result = ALLOC_ITEM(ArrayList_s, "ArrayList structure");
    if (Result != NULL)
    {
        Result->Array           = NULL;
        Result->Type            = Type;
        Result->ItemSize        = GetElementSize(Type);
        Result->Count           = 0;
        Result->Capacity        = 0;
        Result->IsVisitingItems = FALSE;
        if (CapacityRequestAdjuster != NULL)
        {
            /* install the client's capacity request adjuster */
            Result->CapacityRequestAdjuster           = CapacityRequestAdjuster;
            Result->CapacityRequestAdjusterClientData = CapacityRequestAdjusterClientData;
        }
        else
        {
            /* install the default capacity request adjuster */
            Result->CapacityRequestAdjuster           = AdjustCapacityRequest;
            Result->CapacityRequestAdjusterClientData = 0;
        }

        /* enalarge the list to the estimated capacity */
        if (!ArrayListEnlargeCapacity(Result, EstimatedCapacity))
            ArrayListDealloc(&Result, NULL, 0);
    }

    ENSURE(ArrayListIsValid(Result) || Result == NULL);
    ENSURE(IMPLICATION(Result != NULL, Result->Capacity >= EstimatedCapacity));
    return Result;
}


/**
 * Deallocates the list handle and set the handle to NULL.
 *
 * param ArrayList
 *     Reference to an array list handle.
 * param ItemDestructor
 *     Destructor responsible for array list item cleanup or
 *     NULL if no item cleanup is desired.
 * param ClientData
 *     Any client data needed for cleanup.
 */
void ArrayListDealloc(ArrayList_pa               *ArrayList,
                      ArrayListItemDestructor_pf  ItemDestructor,
                      ArbParam_t                  ClientData)
{
    REQUIRE(VALID_REF(ArrayList));
    REQUIRE(ArrayListIsValid(*ArrayList) || *ArrayList == NULL);
    REQUIRE(VALID_FN_REF(ItemDestructor) || ItemDestructor == NULL);

    if (*ArrayList != NULL)
    {
        /* request item cleanup if a destructor was supplied */
        if (ItemDestructor != NULL && (*ArrayList)->Count != 0)
            DestroyItems(*ArrayList, 0, (*ArrayList)->ItemSize, (*ArrayList)->Count,
                         ItemDestructor, ClientData);

        /* release the list */
        if ((*ArrayList)->Capacity != 0)
            FREE_ARRAY((*ArrayList)->Array, "array list");

        /* release the list structure itself */
        FREE_ITEM(*ArrayList, "ArrayList structure");
        *ArrayList = NULL;
    }

    ENSURE(*ArrayList == NULL);
}


#if !defined USE_MACROS_FOR_FUNCTIONS
/**
 * Gets the number of items currently maintained by the list.
 *
 * param
 *     Array list in question.
 *
 * return
 *     Number of items maintained by the list.
 */
LgIndex_t ArrayListGetCount_FUNC(ArrayList_pa ArrayList)
{
    LgIndex_t Result;

    REQUIRE(ArrayListIsValid(ArrayList));

    Result = ArrayListGetCount_MACRO(ArrayList);

    ENSURE(Result >= 0);
    return Result;
}
#endif


/**
 * Empties the array list of all items.
 *
 * param ArrayList
 *     Array list from which to delete all items.
 * param ItemDestructor
 *     Destructor responsible for array list item cleanup or
 *     NULL if no item cleanup is desired.
 * param ClientData
 *     Any client data needed for cleanup.
 */
void ArrayListDeleteAllItems(ArrayList_pa               ArrayList,
                             ArrayListItemDestructor_pf ItemDestructor,
                             ArbParam_t                 ClientData)
{
    REQUIRE(ArrayListIsValid(ArrayList));
    REQUIRE(VALID_FN_REF(ItemDestructor) || ItemDestructor == NULL);
    REQUIRE(!ArrayList->IsVisitingItems);

    /* request item cleanup if a destructor was supplied */
    if (ItemDestructor != NULL && ArrayList->Count != 0)
        DestroyItems(ArrayList, 0, ArrayList->ItemSize, ArrayList->Count,
                     ItemDestructor, ClientData);

    /*
     * Fill the vacated items with zeros. This default value of zero is necessary
     * for many other array list operations.
     */
    if (ArrayList->Count != 0)
        memset(&ArrayList->Array[0], 0, ArrayList->Count*ArrayList->ItemSize);

    ArrayList->Count = 0;

    ENSURE(ArrayListIsValid(ArrayList) && ArrayList->Count == 0);
}


/**
 * Deletes 'Count' items from the array list. The members following the
 * items deleted are shifted down accordingly to fill the vacated space.
 *
 * param ArrayList
 *     Array list containing the items to delete.
 * param ItemOffset
 *     Offset to the first item to delete in the list.
 * param Count
 *     Number of items to delete.
 * param ItemDestructor
 *     Destructor responsible for array list item cleanup or
 *     NULL if no item cleanup is desired.
 * param ClientData
 *     Any client data needed for cleanup.
 */
void ArrayListDeleteItems(ArrayList_pa               ArrayList,
                          LgIndex_t                  ItemOffset,
                          LgIndex_t                  Count,
                          ArrayListItemDestructor_pf ItemDestructor,
                          ArbParam_t                 ClientData)
{
    REQUIRE(ArrayListIsValid(ArrayList));
    REQUIRE(0 <= ItemOffset && ItemOffset <= ArrayList->Count - 1);
    REQUIRE(1 <= Count && ItemOffset + Count <= ArrayList->Count);
    REQUIRE(VALID_FN_REF(ItemDestructor) || ItemDestructor == NULL);
    REQUIRE(!ArrayList->IsVisitingItems);

    /* release the items if a destructor is installed */
    if (ItemDestructor != NULL)
        DestroyItems(ArrayList, ItemOffset, ArrayList->ItemSize, Count,
                     ItemDestructor, ClientData);

    /* if we deleted the items from the middle of the array then     */
    /* shift the end items down by 'Count' to fill the vacated space */
    if (ItemOffset + Count < ArrayList->Count)
        CopyArrayItems(ArrayList->Array, ItemOffset,
                       ArrayList->Array, ItemOffset + Count,
                       ArrayList->Count - (ItemOffset + Count),
                       ArrayList->ItemSize);
    /*
     * Fill the vacated items with zeros. This default value of zero is necessary
     * for many other array list operations.
     */
    memset(&ArrayList->Array[(ArrayList->Count - Count)*ArrayList->ItemSize],
           0, Count*ArrayList->ItemSize);

    /* update the count but leave the capacity alone */
    ArrayList->Count -= Count;

    ENSURE(ArrayListIsValid(ArrayList));
}


/**
 * Deletes an item from the array list. The members following the item
 * deleted are shifted down accordingly to fill the vacated space.
 *
 * param ArrayList
 *     Array list containing the item to delete.
 * param ItemOffset
 *     Offset to the item in the list.
 * param ItemDestructor
 *     Destructor responsible for array list item cleanup or
 *     NULL if no item cleanup is desired.
 * param ClientData
 *     Any client data needed for cleanup.
 */
void ArrayListDeleteItem(ArrayList_pa               ArrayList,
                         LgIndex_t                  ItemOffset,
                         ArrayListItemDestructor_pf ItemDestructor,
                         ArbParam_t                 ClientData)
{
    REQUIRE(ArrayListIsValid(ArrayList));
    REQUIRE(0 <= ItemOffset && ItemOffset <= ArrayList->Count - 1);
    REQUIRE(VALID_FN_REF(ItemDestructor) || ItemDestructor == NULL);

    ArrayListDeleteItems(ArrayList, ItemOffset, 1, ItemDestructor, ClientData);

    ENSURE(ArrayListIsValid(ArrayList));
}


/**
 * Removes 'Count' items from the array list beginning at the specified
 * item offset. The members following the items removed are shifted down
 * accordingly to fill the vacated space.
 *
 * param ArrayList
 *     Array list containing the items to remove.
 * param ItemOffset
 *     Offset to the first item to remove in the list.
 * param Count
 *     Number of items to remove.
 *
 * return
 *     Array list handle referring to the removed items if sufficient
 *     memory was available, otherwise a handle to NULL.
 */
ArrayList_pa ArrayListRemoveItems(ArrayList_pa ArrayList,
                                  LgIndex_t    ItemOffset,
                                  LgIndex_t    Count)
{
    ArrayList_pa Result;

    REQUIRE(ArrayListIsValid(ArrayList));
    REQUIRE(0 <= ItemOffset && ItemOffset <= ArrayList->Count - 1);
    REQUIRE(1 <= Count && ItemOffset + Count <= ArrayList->Count);
    REQUIRE(!ArrayList->IsVisitingItems);

    /* get a copy of the items and delete them from the source */
    Result = ArrayListGetItems(ArrayList, ItemOffset, Count);
    if (Result != NULL)
        ArrayListDeleteItems(ArrayList, ItemOffset, Count, NULL, 0);

    ENSURE(ArrayListIsValid(ArrayList));
    ENSURE(ArrayListIsValid(Result) || Result == NULL);
    return Result;
}


/**
 * Removes an item from the array list. The members following the item
 * removed are shifted down accordingly to fill the vacated space.
 *
 * param ArrayList
 *     Array list containing the item to remove.
 * param ItemOffset
 *     Offset to the item in the list.
 *
 * return
 *     Item removed from the array list.
 */
ArrayListItem_u ArrayListRemoveItem(ArrayList_pa ArrayList,
                                    LgIndex_t    ItemOffset)
{
    ArrayListItem_u Result;

    REQUIRE(ArrayListIsValid(ArrayList));
    REQUIRE(0 <= ItemOffset && ItemOffset <= ArrayList->Count - 1);
    REQUIRE(!ArrayList->IsVisitingItems);

    /* record the original item */
    CopyArrayItems((char *)&Result, 0,
                   ArrayList->Array, ItemOffset,
                   1, ArrayList->ItemSize);

    /* delete the item from the array */
    ArrayListDeleteItems(ArrayList, ItemOffset, 1, NULL, 0);

    ENSURE(ArrayListIsValid(ArrayList));
    return Result;
}


/**
 * Inserts copies of the items from the source list to the target list at
 * the specified offset. The target list will expand to accommodate the
 * additional items. The source list remains unchanged.
 *
 * param Target
 *     Array list receiving the source items.
 * param ItemOffset
 *     Offset at which to insert the source list items.
 * param Source
 *     Array list supplying the source items.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t ArrayListInsert(ArrayList_pa Target,
                          LgIndex_t    ItemOffset,
                          ArrayList_pa Source)
{
    Boolean_t IsOk = TRUE;

    REQUIRE(ArrayListIsValid(Target));
    REQUIRE(ItemOffset >= 0);
    REQUIRE(ArrayListIsValid(Source));
    REQUIRE(Target != Source);
    REQUIRE(Target->Type == Source->Type);
    REQUIRE(!Target->IsVisitingItems);

    if (Source->Count != 0)
    {
        LgIndex_t NeededCapacity;

        /* if necessary enlarge the target list to accommodate the request */
        if (ItemOffset > Target->Count)
            NeededCapacity = ItemOffset + Source->Count;
        else
            NeededCapacity = Target->Count + Source->Count;
        if (NeededCapacity > Target->Capacity)
            IsOk = ArrayListEnlargeCapacity(Target, NeededCapacity);

        if (IsOk)
        {
            if (ItemOffset < Target->Count)
            {
                /* shift all items in the target list ahead of the  */
                /* insert position up by the number of items in the */
                /* source list to make room for the new items       */
                CopyArrayItems(Target->Array, ItemOffset + Source->Count,
                               Target->Array, ItemOffset,
                               Target->Count - ItemOffset,
                               Target->ItemSize);
                Target->Count += Source->Count;
            }
            else
            {
                /* no shifting to do, just update the count */
                if (ItemOffset > Target->Count)
                    Target->Count = ItemOffset + Source->Count;
                else
                    Target->Count += Source->Count;
            }

            /* insert the items and update the count */
            CopyArrayItems(Target->Array, ItemOffset,
                           Source->Array, 0,
                           Source->Count, Source->ItemSize);
        }
    }

    ENSURE(ArrayListIsValid(Target));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}


/**
 * Inserts the item into the array list at the specified offset. The list will
 * be expanded to accommodate the additional item. If the offset is beyond the
 * end of the list it is sized accordingly and the intervening items between
 * the last item of the original state and the last item of the new state are
 * guaranteed to be 0.
 *
 * param ArrayList
 *     Array list target in which to insert the item.
 * param ItemOffset
 *     Offset at which to insert the item.
 * param Item
 *     Item to insert at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t ArrayListInsertItem(ArrayList_pa    ArrayList,
                              LgIndex_t       ItemOffset,
                              ArrayListItem_u Item)
{
    Boolean_t IsOk = TRUE;
    LgIndex_t NeededCapacity;

    REQUIRE(ArrayListIsValid(ArrayList));
    REQUIRE(ItemOffset >= 0);
    REQUIRE(!ArrayList->IsVisitingItems);

    /* if necessary enlarge the list to accommodate the request */
    if (ItemOffset > ArrayList->Count)
        NeededCapacity = ItemOffset + 1;
    else
        NeededCapacity = ArrayList->Count + 1;
    if (NeededCapacity > ArrayList->Capacity)
        IsOk = ArrayListEnlargeCapacity(ArrayList, NeededCapacity);

    if (IsOk)
    {
        if (ItemOffset < ArrayList->Count)
        {
            /* shift all items in the target list ahead of the insert */
            /* position up by one to make room for the new item       */
            CopyArrayItems(ArrayList->Array, ItemOffset + 1,
                           ArrayList->Array, ItemOffset,
                           ArrayList->Count - ItemOffset,
                           ArrayList->ItemSize);
            ArrayList->Count++;
        }
        else
        {
            /* no shifting to do, just update the count */
            if (ItemOffset > ArrayList->Count)
                ArrayList->Count = ItemOffset + 1;
            else
                ArrayList->Count++;
        }

        /* insert the item */
        CopyArrayItems(ArrayList->Array, ItemOffset,
                       (char *)&Item, 0,
                       1, ArrayList->ItemSize);
    }

    ENSURE(ArrayListIsValid(ArrayList));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}


/**
 * Visits array list items calling the item visitor for each item.
 *
 * param ArrayList
 *     Array list needing its items destroyed.
 * param ItemOffset
 *     Offset to the first item to visit from the list.
 * param Count
 *     Number of items to visit.
 * param ItemVisitor
 *     Function called for each array list item.
 * param CientData
 *     Any client data needed for the visitor.
 *
 * return
 *     TRUE if the all element were visited, otherwise
 *     FALSE if the visitation was terminated early
 */
Boolean_t ArrayListVisitItems(ArrayList_pa            ArrayList,
                              LgIndex_t               ItemOffset,
                              LgIndex_t               Count,
                              ArrayListItemVisitor_pf ItemVisitor,
                              ArbParam_t              ClientData)
{
    Boolean_t   DoContinue = TRUE;
    Boolean_t   IsVisitingItems;
    SmInteger_t ItemSize;
    LgIndex_t   Index;

    REQUIRE(ArrayListIsValid(ArrayList));
    REQUIRE(VALID_FN_REF(ItemVisitor));

    IsVisitingItems = ArrayList->IsVisitingItems;
    ArrayList->IsVisitingItems = TRUE; /* guards against structure changes */

    for (Index = 0, ItemSize = ArrayList->ItemSize;
         Index < Count && DoContinue;
         Index++)
    {
        LgIndex_t ItemIndex = (Index + ItemOffset) * ItemSize;
        DoContinue = ItemVisitor((void *) & ArrayList->Array[ItemIndex], ClientData);
    }

    ArrayList->IsVisitingItems = IsVisitingItems;

    ENSURE(ArrayList->IsVisitingItems == IsVisitingItems);
    ENSURE(VALID_BOOLEAN(DoContinue));
    return DoContinue;
}


/**
 * Gets copies of 'Count' items from the array list beginning at the
 * specified item offset. Note that if the items are pointer types
 * the copies are of the pointers and not the pointees.
 *
 * param ArrayList
 *     Array list containing the items to copy.
 * param ItemOffset
 *     Offset to the first item to copy from the list.
 * param Count
 *     Number of items to copy.
 *
 * return
 *     Array list handle referring to the copied items if sufficient
 *     memory was available, otherwise a handle to NULL.
 */
ArrayList_pa ArrayListGetItems(ArrayList_pa ArrayList,
                               LgIndex_t    ItemOffset,
                               LgIndex_t    Count)
{
    ArrayList_pa Result;

    REQUIRE(ArrayListIsValid(ArrayList));
    REQUIRE(0 <= ItemOffset && ItemOffset <= ArrayList->Count - 1);
    REQUIRE(1 <= Count && ItemOffset + Count <= ArrayList->Count);

    Result = ArrayListAlloc(Count, ArrayList->Type,
                            ArrayList->CapacityRequestAdjuster,
                            ArrayList->CapacityRequestAdjusterClientData);
    if (Result != NULL)
    {
        /* copy the original items into the result */
        CopyArrayItems(Result->Array, 0,
                       ArrayList->Array, ItemOffset,
                       Count, ArrayList->ItemSize);
        Result->Count = Count;
    }

    ENSURE(ArrayListIsValid(ArrayList));
    ENSURE(ArrayListIsValid(Result) || Result == NULL);
    return Result;
}


/**
 * Gets the item at the specified offset in the list.
 *
 * param ArrayList
 *     Array list containing the desired item.
 * param ItemOffset
 *     Offset to the item in the list.
 *
 * return
 *     The requested item.
 */
ArrayListItem_u ArrayListGetItem(ArrayList_pa ArrayList,
                                 LgIndex_t    ItemOffset)
{
    ArrayListItem_u Result;

    REQUIRE(ArrayListIsValid(ArrayList));
    REQUIRE(0 <= ItemOffset && ItemOffset <= ArrayList->Count - 1);

    CopyArrayItems((char *)&Result, 0,
                   ArrayList->Array, ItemOffset,
                   1, ArrayList->ItemSize);

    return Result;
}


#if !defined USE_MACROS_FOR_FUNCTIONS
/**
 * Gets the item's internal reference at the specified offset in the list.
 *
 * WARNING:
 *     Some array list functions modify the internal buffer.
 *     This will invalidate this reference however it is
 *     the client's responsibility not to make further use
 *     of it. In addition, this reference should never be
 *     deallocated directly as the array list assumes the
 *     responsible for the cleanup.
 *
 * param ArrayList
 *     Array list containing the desired item.
 * param ItemOffset
 *     Offset to the item in the list.
 *
 * return
 *     The internal reference to the requested item.
 */
const void *ArrayListGetItemInternalRef_FUNC(ArrayList_pa ArrayList,
                                             LgIndex_t    ItemOffset)
{
    REQUIRE(ArrayListIsValid(ArrayList));
    REQUIRE(0 <= ItemOffset && ItemOffset <= ArrayList->Count - 1);

    const void *Result = ArrayListGetItemInternalRef_MACRO(ArrayList, ItemOffset);
    ENSURE(Result == NULL || VALID_REF(Result));
    return Result;
}
#endif


/**
 * Places the item at the specified offset. If the offset is beyond the
 * end of the list it is sized accordingly and the intervening items
 * between the last item of the original state and the last item of the
 * new state are guaranteed to be 0.
 *
 * param ArrayList
 *     Array list target in which to set the item.
 * param ItemOffset
 *     Offset of the item.
 * param Item
 *     Item to set at the specified offset.
 * param ItemDestructor
 *     Destructor responsible for array list item cleanup or
 *     NULL if no item cleanup is desired.
 * param ClientData
 *     Any client data needed for cleanup.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t ArrayListSetItem(ArrayList_pa               ArrayList,
                           LgIndex_t                  ItemOffset,
                           ArrayListItem_u            Item,
                           ArrayListItemDestructor_pf ItemDestructor,
                           ArbParam_t                 ClientData)
{
    Boolean_t IsOk = TRUE;

    REQUIRE(ArrayListIsValid(ArrayList));
    REQUIRE(ItemOffset >= 0);
    REQUIRE(VALID_FN_REF(ItemDestructor) || ItemDestructor == NULL);
    REQUIRE(IMPLICATION(ItemOffset + 1 > ArrayList->Count,
                        !ArrayList->IsVisitingItems));

    /* release the item if a destructor is installed */
    if (ItemDestructor != NULL && ItemOffset < ArrayList->Count)
        DestroyItems(ArrayList, ItemOffset, ArrayList->ItemSize, 1,
                     ItemDestructor, ClientData);

    /* if necessary enlarge the list to accommodate the request */
    if (ItemOffset + 1 > ArrayList->Capacity)
        IsOk = ArrayListEnlargeCapacity(ArrayList, ItemOffset + 1);

    if (IsOk)
    {
        if (ItemOffset + 1 > ArrayList->Count)
            ArrayList->Count = ItemOffset + 1;
        CopyArrayItems(ArrayList->Array, ItemOffset,
                       (char *)&Item, 0,
                       1, ArrayList->ItemSize);
    }

    ENSURE(ArrayListIsValid(ArrayList));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}


/**
 * Appends the item to the list. The list will be expanded
 * to accommodate the additional item.
 *
 * param ArrayList
 *     Array list target to which the item is to be appended.
 * param Item
 *     Item to append to the array list.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t ArrayListAppendItem(ArrayList_pa    ArrayList,
                              ArrayListItem_u Item)
{
    Boolean_t IsOk;

    REQUIRE(ArrayListIsValid(ArrayList));
    REQUIRE(!ArrayList->IsVisitingItems);

    IsOk = ArrayListInsertItem(ArrayList, ArrayList->Count, Item);

    ENSURE(ArrayListIsValid(ArrayList));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}


/**
 * Appends copies of the items from the source list to the target list.
 * The source list remains unchanged.
 *
 * param Target
 *     Array list receiving the source items.
 * param Source
 *     Array list supplying the source items.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t ArrayListAppend(ArrayList_pa Target,
                          ArrayList_pa Source)
{
    Boolean_t IsOk;

    REQUIRE(ArrayListIsValid(Target));
    REQUIRE(ArrayListIsValid(Source));
    REQUIRE(Target != Source);
    REQUIRE(Target->Type == Source->Type);
    REQUIRE(!Target->IsVisitingItems);

    IsOk = ArrayListInsert(Target, Target->Count, Source);

    ENSURE(ArrayListIsValid(Target));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}


/**
 * Copies the items of the array list.
 *
 * param ArrayList
 *     Array list to copy.
 * param ItemDuplicator
 *     Duplicator responsible for array list item duplication or
 *     NULL if an exact item copy is desired. In other words for
 *     pointer types the effect is copy by reference.
 * param ClientData
 *     Any client data needed for duplication.
 *
 * return
 *     Handle to a duplicate of the specified array list if sufficient
 *     memory permitted the operation, otherwise NULL.
 */
ArrayList_pa ArrayListCopy(ArrayList_pa               ArrayList,
                           ArrayListItemDuplicator_pf ItemDuplicator,
                           ArbParam_t                 ClientData)
{
    ArrayList_pa Result;

    REQUIRE(ArrayListIsValid(ArrayList));
    REQUIRE(VALID_FN_REF(ItemDuplicator) || ItemDuplicator == NULL);

    Result = ArrayListAlloc(ArrayList->Count, ArrayList->Type,
                            ArrayList->CapacityRequestAdjuster,
                            ArrayList->CapacityRequestAdjusterClientData);
    if (Result != NULL && ArrayList->Count != 0)
    {
        Boolean_t IsOk = TRUE;
        if (ItemDuplicator != NULL)
            /* client defines how the item duplication is performed */
            IsOk = DuplicateItems(Result->Array, 0,
                                  ArrayList->Array, 0,
                                  ArrayList->ItemSize, ArrayList->Count,
                                  ItemDuplicator, ClientData);
        else
            /* copy the original items into the result */
            CopyArrayItems(Result->Array, 0,
                           ArrayList->Array, 0,
                           ArrayList->Count,
                           ArrayList->ItemSize);
        if (IsOk)
            Result->Count = ArrayList->Count;
        else
            ArrayListDealloc(&Result, NULL, 0);
    }

    ENSURE(Result == NULL ||
           (ArrayListIsValid(Result) && Result->Count == ArrayList->Count));
    return Result;
}


/**
 * Creates a native 'C' array containing copies of the items held in the
 * source array list.
 *
 * param Source
 *     Array list containing the items of interest.
 * param ItemDuplicator
 *     Duplicator responsible for array list item duplication or
 *     NULL if an exact item copy is desired. In other words for
 *     pointer types the effect is copy by reference.
 * param ClientData
 *     Any client data needed for duplication.
 *
 * return
 *     Allocated array populated with copies of the members of the array list
 *     or NULL if there are no items in the list or if the allocation was
 *     not successful. The caller is responsible for deallocation of the
 *     array (but not the individual members unless a item duplication function
 *     was supplied) when it is no longer needed.
 */
void *ArrayListToArray(ArrayList_pa               Source,
                       ArrayListItemDuplicator_pf ItemDuplicator,
                       ArbParam_t                 ClientData)
{
    void *Result;

    REQUIRE(ArrayListIsValid(Source));
    REQUIRE(VALID_FN_REF(ItemDuplicator) || ItemDuplicator == NULL);

    if (Source->Count != 0)
        Result = (void *)ALLOC_ARRAY(Source->Count * Source->ItemSize,
                                     char, "native array");
    else
        Result = NULL;

    if (Result != NULL)
    {
        Boolean_t IsOk = TRUE;
        if (ItemDuplicator != NULL)
            /* client defines how the item duplication is performed */
            IsOk = DuplicateItems((char*)Result, 0,
                                  Source->Array, 0,
                                  Source->ItemSize, Source->Count,
                                  ItemDuplicator, ClientData);
        else
            /* copy the original items into the result */
            CopyArrayItems((char *)Result, 0,
                           Source->Array, 0,
                           Source->Count,
                           Source->ItemSize);
        if (!IsOk)
        {
            /* Hack to remove delete warning... */
            char *Tmp = (char *)Result;
            FREE_ARRAY(Tmp, "native array");
        }
    }

    ENSURE(VALID_REF(Result) || Result == NULL);
    return Result;
}


/**
 * Creates an array list containing copies of the items held in the
 * native 'C' array.
 *
 * param Source
 *     Native 'C' array containing the items of interest.
 * param Count
 *     Number of items contained in the native 'C' array.
 * param Type
 *     Type of items contained in the native 'C' array.
 * param ItemDuplicator
 *     Duplicator responsible for array list item duplication or
 *     NULL if an exact item copy is desired. In other words for
 *     pointer types the effect is copy by reference.
 * param ClientData
 *     Any client data needed for duplication.
 *
 * return
 *     Array list handle containing copies of the items held in the
 *     native 'C' array if sufficient memory was available, otherwise
 *     a handle to NULL.
 */
ArrayList_pa ArrayListFromArray(void                       *Source,
                                LgIndex_t                   Count,
                                ArrayListType_e             Type,
                                ArrayListItemDuplicator_pf  ItemDuplicator,
                                ArbParam_t                  ClientData)
{
    ArrayList_pa Result;

    REQUIRE(VALID_REF(Source));
    REQUIRE(Count >= 0);
    REQUIRE(VALID_ENUM(Type, ArrayListType_e));
    REQUIRE(VALID_FN_REF(ItemDuplicator) || ItemDuplicator == NULL);

    Result = ArrayListAlloc(Count, Type, NULL, 0);
    if (Result != NULL && Count != 0)
    {
        Boolean_t IsOk = TRUE;
        if (ItemDuplicator != NULL)
            /* client defines how the item duplication is performed */
            IsOk = DuplicateItems(Result->Array, 0,
                                  (char*)Source, 0,
                                  Result->ItemSize, Count,
                                  ItemDuplicator, ClientData);
        else
            /* copy the original items into the result */
            CopyArrayItems(Result->Array, 0,
                           (char *)Source, 0,
                           Count, Result->ItemSize);
        if (IsOk)
            Result->Count = Count;
        else
            ArrayListDealloc(&Result, NULL, 0);
    }

    ENSURE(ArrayListIsValid(Result) || Result == NULL);
    return Result;
}


/**
 * Holds the comparator function pointer for sorting.
 */
static ArrayListItemComparator_pf ComparatorFunction = NULL;


/**
 * Holds the context for comparisons. This information is forwarded to
 * the item comparator function for sorting.
 */
static ArbParam_t ComparatorClientData = 0;


/**
 * Holds the item size of the individual array components for sorting.
 */
static SmInteger_t ComparatorItemSize = 0;


/**
 * Forwards the comparison test to the 'Comparator' supplied to the
 * 'ArrayListQSort' function.
 *
 * param Item1Ref
 *     Reference to base address of Item1.
 * param Item2Ref
 *     Reference to base address of Item2.
 *
 * return
 *     - A value less than zero if Item1 is less than Item2.
 *     - A value of zero if Item1 is equal to Item2.
 *     - A value greater than zero if Item1 is greater than Item2.
 */
static int QSortCompareItemsAdapter(const void *Item1Ref,
                                    const void *Item2Ref)
{
    int Result;
    ArrayListItem_u Item1;
    ArrayListItem_u Item2;

    REQUIRE(Item1Ref != NULL);
    REQUIRE(Item2Ref != NULL);

    /* collect up the items */
    CopyArrayItems((char *)&Item1, 0,
                   (char *)Item1Ref, 0,
                   1, ComparatorItemSize);
    CopyArrayItems((char *)&Item2, 0,
                   (char *)Item2Ref, 0,
                   1, ComparatorItemSize);

    /* forward the call */
    Result = ComparatorFunction(Item1, Item2, ComparatorClientData);

    ENSURE(Result == -1 || Result == 0 || Result == 1);
    return Result;
}

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

/**
 * Sorts the array list using the qsort algorithm.
 *
 * param ArrayList
 *     Array list to sort.
 * param Comparator
 *     Function called to compare two array list elements.
 * param ClientData
 *     Contextual information that is passed along to the comparator function.
 */
void ArrayListQSort(ArrayList_pa                ArrayList,
                    ArrayListItemComparator_pf  Comparator,
                    ArbParam_t                  ClientData)
{
    ArrayListItemComparator_pf CurComparatorFunction;
    ArbParam_t                 CurComparatorClientData;
    SmInteger_t                CurComparatorItemSize;

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    REQUIRE(ArrayListIsValid(ArrayList));
    REQUIRE(VALID_FN_REF(Comparator));

    /* to support sort recursion we need to save off the current values */
    CurComparatorFunction   = ComparatorFunction;
    CurComparatorClientData = ComparatorClientData;
    CurComparatorItemSize   = ComparatorItemSize;

    /* set up for comparison proxy */
    ComparatorFunction   = Comparator;
    ComparatorClientData = ClientData;
    ComparatorItemSize   = ArrayList->ItemSize;

    /* sort the array */
    qsort(ArrayList->Array, ArrayList->Count,
          ArrayList->ItemSize, QSortCompareItemsAdapter);

    /* cleanup */
    ComparatorFunction   = CurComparatorFunction;
    ComparatorClientData = CurComparatorClientData;
    ComparatorItemSize   = CurComparatorItemSize;

    ENSURE(ArrayListIsValid(ArrayList));
}

/**
 * Binary searches a sorted array looking for a match using the supplied
 * comparator function. If a match is found the resulting item index refers to
 * the location otherwise it refers to the location where the item could be
 * inserted in sorted order.
 *
 * param ArrayList
 *     Array list to sort.
 * param Item
 *     The item for which to search.
 * param Comparator
 *     Function called to compare the Item to the array list elements.
 * param ClientData
 *     Contextual information that is passed along to the comparator function.
 * param ItemIndex
 *     Pointer to the resulting position where the item was found or where the
 *     item could be inserted in sorted order. If the pointer is NULL the
 *     result position is not returned.
 *
 * result
 *     TRUE if the item was found in the list, FALSE otherwise.
 */
Boolean_t ArrayListBSearch(ArrayList_pa                ArrayList,
                           ArrayListItem_u             Item,
                           ArrayListItemComparator_pf  Comparator,
                           ArbParam_t                  ClientData,
                           LgIndex_t                  *ItemIndex)
{
    REQUIRE(ArrayListIsValid(ArrayList));
    REQUIRE(VALID_FN_REF(Comparator));
    REQUIRE(ItemIndex == NULL || VALID_REF(ItemIndex));

    LgIndex_t MiddleItemIndex = 0;
    LgIndex_t FirstItemIndex  = 0;
    LgIndex_t NumItems        = ArrayListGetCount(ArrayList);
    LgIndex_t LastItemIndex   = NumItems - 1;
    Boolean_t Found = FALSE;
    while (FirstItemIndex <= LastItemIndex && !Found)
    {
        /* calculate the middle item index for current sub-range */
        MiddleItemIndex = (FirstItemIndex + LastItemIndex) / 2;

        int CompareResult = Comparator(ArrayListGetItem(ArrayList, MiddleItemIndex), Item, ClientData);
        if (CompareResult > 0)
            LastItemIndex = MiddleItemIndex - 1;
        else if (CompareResult < 0)
            FirstItemIndex = MiddleItemIndex + 1;
        else
            Found = TRUE;
    }

    if (ItemIndex != NULL)
    {
        if (Found || NumItems == 0 || FirstItemIndex < NumItems)
            *ItemIndex = MiddleItemIndex;
        else
            *ItemIndex = NumItems; /* ...in other words it goes on the end */
    }

    ENSURE(IMPLICATION(ItemIndex != NULL,
                       0 <= *ItemIndex && *ItemIndex <= ArrayListGetCount(ArrayList)));
    ENSURE(VALID_BOOLEAN(Found));
    return Found;
}

#if !defined USE_MACROS_FOR_FUNCTIONS
/**
 * Gets the array list's internal buffer representation.
 *
 * WARNING:
 *     Some array list functions modify the internal buffer.
 *     This will invalidate this reference however it is
 *     the client's responsibility not to make further use
 *     of it. In addition, this reference should never be
 *     deallocated directly as the array list assumes the
 *     responsible for the cleanup.
 *
 * param ArrayList
 *     Array list for which a reference to the internal
 *     buffer is desired.
 *
 * return
 *     Reference to the array list's internal buffer.
 */
const void *ArrayListGetInternalRef_FUNC(ArrayList_pa ArrayList)
{
    REQUIRE(ArrayListIsValid(ArrayList));

    const void *Result = ArrayListGetInternalRef_MACRO(ArrayList);
    ENSURE(Result == NULL || VALID_REF(Result));
    return Result;
}
#endif
